# Regional-inequalities-in-service-provision-and-Risperidone-use-for-children-and-adolescents-with-ASD
################################################################################
# TÍTULO: Análise de Tendência do Consumo de Risperidona (DDD) no Brasil
# AUTOR(A): Alice Ramos (adaptado por Gemini)
# DATA: 30 de outubro de 2025
#
# DESCRIÇÃO:
# Este script analisa os dados de dispensação de risperidona para pacientes
# com TEA (0-19 anos), calcula a Dose Diária Definida (DDD), normaliza por
# população e realiza uma análise de tendência temporal (Regressão de 
# Prais-Winsten) por região do Brasil.
#
################################################################################


# === 0. CONFIGURAÇÃO INICIAL ===

# --- 0.1 Limpeza do Ambiente ---
# Limpa todas as variáveis, libera memória e fecha gráficos abertos
rm(list = ls())
gc()
if (!is.null(dev.list())) dev.off()

# --- 0.2 Carregamento de Pacotes ---
# Instala e carrega as bibliotecas necessárias.
# O 'pacman' gerencia a instalação automaticamente.
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  dplyr,      # Manipulação de dados
  stringr,    # Manipulação de texto (strings)
  ggplot2,    # Visualização de dados
  lubridate,  # Manipulação de datas
  tibble,     # Criação de data frames (tribble)
  prais,      # Regressão de Prais-Winsten
  nlme,       # Modelos lineares (fallback para Prais)
  Kendall     # Teste de Mann-Kendall (para exploração)
)

# --- 0.3 Diretório de Trabalho ---
# (!!) RECOMENDAÇÃO PARA O GITHUB (!!)
# É uma má prática usar setwd() em scripts compartilháveis.
# A melhor abordagem é usar "RStudio Projects" (.Rproj).
# O script irá rodar a partir da raiz do projeto.
#
# setwd("D:/") # <--- Linha original (comentada)
#
# Seus arquivos devem estar no diretório do projeto ou em subpastas, por exemplo:
# "dados/dados risperidona TEA.csv"
# "dados/pop 0 19 anos regiao.csv"
# "dados/dados DDD trimestre.csv"

# --- 0.4 Constantes Globais ---
# Definir a DDD da Risperidona em gramas (5 mg = 0.005 g)
DDD_g <- 0.005


# === 1. CARGA E LIMPEZA DOS DADOS ===

# Carrega os dados de dispensação.
# A codificação "Latin1" é usada para tratar caracteres especiais.
dados <- read.csv("dados risperidona TEA.csv",
                  fileEncoding = "Latin1")

# (Opcional) Conferir nomes das colunas e medicamentos únicos
# unique(dados$nome_medicamento)
# colnames(dados)


# === 2. PADRONIZAÇÃO DE DOSES (COMPRIMIDOS VS SOLUÇÃO) ===

# Esta seção extrai a dosagem (em gramas) de diferentes apresentações
# (comprimidos "mg" vs. solução "mg/mL").

dados <- dados %>%
  mutate(
    # Detecta se é solução oral (pela unidade "mg/mL")
    eh_solucao = str_detect(nome_medicamento, "mg/mL"),

    # Para COMPRIMIDOS: extrai "mg" e converte para "g"
    mg_valor = ifelse(!eh_solucao,
                      as.numeric(str_extract(nome_medicamento, "\\d+\\.*\\d*")),
                      NA),
    dose_gramas = mg_valor / 1000, # mg -> g

    # Para SOLUÇÃO ORAL: extrai "mg/mL" e volume do frasco "mL"
    mg_por_ml = ifelse(eh_solucao,
                       as.numeric(str_extract(nome_medicamento, "\\d+")), # Extrai o "1" de "1 mg/mL"
                       NA),
    volume_ml = ifelse(eh_solucao,
                       as.numeric(str_extract(nome_medicamento, "(?<=frasco )\\d+")), # Extrai o "30" de "frasco 30 mL"
                       NA),

    # Converte a dose da solução para gramas
    g_por_ml = mg_por_ml / 1000, # (mg/mL) -> (g/mL)
    g_total_fr = (mg_por_ml * volume_ml) / 1000 # (g) total no frasco
  )


# === 3. CÁLCULO DE DDD (DEFINED DAILY DOSES) ===

# --- 3.1 Tabela de Referência (Preço e Gramas) ---
# Cria uma tabela "de-para" para normalizar as unidades dispensadas.
# AP_PRIPAL_norm é o código do "Princípio Ativo + Apresentação".
ref_risp <- tibble::tribble(
  ~AP_PRIPAL_norm, ~descricao, ~preco_unit, ~gramas_unit,
  604510012, "Risperidona 1 mg (por comprimido)", 0.10, 0.001, # 1mg
  604510020, "Risperidona 2 mg (por comprimido)", 0.11, 0.002, # 2mg
  604510039, "Risperidona 3 mg (por comprimido)", 0.17, 0.003, # 3mg
  604510047, "Risperidona 1 mg/mL sol. oral (frasco 30 mL)", 21.41, 0.030  # 1mg/mL * 30mL = 30mg
)

# --- 3.2 Cálculo de DDD por Linha ---
# Junta a tabela de referência aos dados e calcula o DDD.
dados_calc <- dados %>%
  # Junta preço e gramagem por unidade (comprimido ou frasco)
  left_join(ref_risp, by = "AP_PRIPAL_norm") %>%
  
  # Calcula a quantidade de unidades dispensadas
  mutate(
    # (Valor Total Gasto) / (Preço Unitário) = Qtd. Dispensada
    qtd_disp = AP_VL_AP / preco_unit,
    
    # Arredonda para o inteiro mais próximo (ex: 1.00001 -> 1)
    qtd_disp_int = as.integer(round(qtd_disp)),
    
    # Gramas totais na linha = (Qtd. Unidades) * (Gramas por Unidade)
    g_total_linha = qtd_disp_int * gramas_unit,
    
    # DDDs na linha = (Gramas Totais) / (Gramas por DDD)
    ddd_linha = g_total_linha / DDD_g
  )


# === 4. ENGENHARIA DE ATRIBUTOS (GEOGRAFIA E TEMPO) ===

# --- 4.1 Criação de Trimestres ---
dados_calc <- dados_calc %>%
  mutate(
    trimestre = case_when(
      mes %in% 1:3 ~ 1,
      mes %in% 4:6 ~ 2,
      mes %in% 7:9 ~ 3,
      mes %in% 10:12 ~ 4
    )
  )

# --- 4.2 Dicionário de UF -> Região ---
uf_regiao <- tibble::tibble(
  estado = c("AC","AL","AP","AM","BA","CE","DF","ES","GO","MA","MT","MS",
             "MG","PA","PB","PR","PE","PI","RJ","RN","RS","RO","RR","SC",
             "SP","SE","TO"),
  regiao = c("Norte","Nordeste","Norte","Norte","Nordeste","Nordeste","Centro-Oeste",
             "Sudeste","Centro-Oeste","Nordeste","Centro-Oeste","Centro-Oeste",
             "Sudeste","Norte","Nordeste","Sul","Nordeste","Nordeste",
             "Sudeste","Nordeste","Sul","Norte","Norte","Sul","Sudeste",
             "Nordeste","Norte")
)

# --- 4.3 Mapeamento Geográfico (Município -> Estado -> Região) ---
# Os 2 primeiros dígitos do código IBGE do município (AP_MUNPCN) são o código do estado.
dados_completo <- dados_calc %>%
  mutate(estado_cod = substr(AP_MUNPCN, 1, 2)) %>%
  
  # Converte o código do estado (IBGE) para a sigla (UF)
  mutate(estado = case_when(
    estado_cod == "12" ~ "AC", estado_cod == "27" ~ "AL", estado_cod == "16" ~ "AP",
    estado_cod == "13" ~ "AM", estado_cod == "29" ~ "BA", estado_cod == "23" ~ "CE",
    estado_cod == "53" ~ "DF", estado_cod == "32" ~ "ES", estado_cod == "52" ~ "GO",
    estado_cod == "21" ~ "MA", estado_cod == "51" ~ "MT", estado_cod == "50" ~ "MS",
    estado_cod == "31" ~ "MG", estado_cod == "15" ~ "PA", estado_cod == "25" ~ "PB",
    estado_cod == "41" ~ "PR", estado_cod == "26" ~ "PE", estado_cod == "22" ~ "PI",
    estado_cod == "33" ~ "RJ", estado_cod == "24" ~ "RN", estado_cod == "43" ~ "RS",
    estado_cod == "11" ~ "RO", estado_cod == "14" ~ "RR", estado_cod == "42" ~ "SC",
    estado_cod == "35" ~ "SP", estado_cod == "28" ~ "SE", estado_cod == "17" ~ "TO",
    TRUE ~ NA_character_
  )) %>%
  
  # Junta a informação da Região
  left_join(uf_regiao, by = "estado")


# === 5. AGREGAÇÃO E NORMALIZAÇÃO POR POPULAÇÃO ===

# --- 5.1 Agregação do DDD por Região/Trimestre ---
# (Opcional: agregado por estado)
# ddd_estado <- dados_completo %>%
#   group_by(estado, trimestre, ano) %>%
#   summarise(ddd_total = sum(ddd_linha, na.rm = TRUE)) %>%
#   arrange(desc(ddd_total))

# Agregado por Região
ddd_regiao <- dados_completo %>%
  group_by(regiao, trimestre, ano) %>%
  summarise(ddd_total = sum(ddd_linha, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(regiao, ano, trimestre)

# --- 5.2 Carga dos Dados Populacionais ---
pop_regiao <- read.csv("pop 0 19 anos regiao.csv", sep=";")

# --- 5.3 Normalização (DDD por 1000 habitantes/dia) ---
dados_final <- ddd_regiao %>%
  left_join(pop_regiao, by = c("regiao", "ano")) %>%
  
  # Define o número de dias em cada trimestre (considerando anos bissextos implicitamente)
  mutate(
    dias_trimestre = case_when(
      trimestre == 1 ~ 90, # Jan-Mar (ou 91 em bissexto, mas 90 é padrão)
      trimestre == 2 ~ 91, # Abr-Jun
      trimestre == 3 ~ 92, # Jul-Set
      trimestre == 4 ~ 92, # Out-Dez
      TRUE ~ 91.25 # Média
    ),
    
    # Fórmula de Normalização:
    # (Total de DDDs no período) / (População * Dias no período) * 1000
    ddd_1000hab_dia = (ddd_total / (Populacao * dias_trimestre)) * 1000
  )

# head(dados_final)


# === 6. UNIFICAÇÃO DOS DADOS (REGIÕES + BRASIL) ===

# Carrega o arquivo separado com os dados já calculados para o "Brasil"
dados_brasil <- read.csv("dados DDD trimestre.csv", sep=";")

# 1) Padroniza colunas das regiões
df_regioes <- dados_final %>%
  select(regiao, ano, trimestre, ddd_1000hab_dia)

# 2) Padroniza colunas do Brasil (adiciona a coluna 'regiao' manualmente)
df_brasil <- dados_brasil %>%
  transmute(
    regiao = "Brasil",
    ano,
    trimestre,
    ddd_1000hab_dia
  )

# 3) Empilha os dataframes (Regiões + Brasil)
df_out <- bind_rows(df_regioes, df_brasil)

# (Opcional) Conferir a estrutura final
# dplyr::glimpse(df_out)


# === 7. VISUALIZAÇÃO DE DADOS (GRÁFICO DE SÉRIE TEMPORAL) ===

# --- 7.1 Preparação para o Gráfico ---
df_out <- df_out %>%
  # Cria uma coluna de data (primeiro dia do trimestre) para o eixo X
  mutate(
    mes_inicio = case_when(
      trimestre == 1 ~ 1,
      trimestre == 2 ~ 4,
      trimestre == 3 ~ 7,
      trimestre == 4 ~ 10
    ),
    data = as.Date(paste(ano, mes_inicio, "01", sep = "-"))
  ) %>%
  
  # Traduz nomes das regiões para inglês (para publicação)
  mutate(regiao = case_when(
    regiao == "Norte" ~ "North",
    regiao == "Nordeste" ~ "Northeast",
    regiao == "Sudeste" ~ "Southeast",
    regiao == "Sul" ~ "South",
    regiao == "Centro-Oeste" ~ "Midwest",
    regiao == "Brasil" ~ "Brazil",
    TRUE ~ regiao
  ))

# --- 7.2 Definição da Paleta de Cores (Acessível) ---
# Paleta Okabe-Ito, amigável para daltonismo
okabe_ito <- c(
  "North"     = "#E69F00", # Laranja
  "Northeast" = "#56B4E9", # Azul claro
  "Southeast" = "#009E73", # Verde
  "South"     = "#F0E442", # Amarelo
  "Midwest"   = "#0072B2", # Azul escuro
  "Brazil"    = "gray10"   # Preto/Cinza
)

# --- 7.3 Plotagem do Gráfico Final ---
ggplot(df_out, aes(x = data, y = ddd_1000hab_dia, color = regiao)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    x = "Year",
    y = "DDD per 1000 inhabitants/day",
    color = "Region"
  ) +
  scale_color_manual(values = okabe_ito) +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y",
    limits = c(as.Date("2016-01-01"), as.Date("2023-12-31"))
  ) +
  scale_y_continuous(
    # Ajusta os "breaks" do eixo Y
    breaks = seq(0, max(df_out$ddd_1000hab_dia, na.rm = TRUE) + 0.025, by = 0.025)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# === 8. ANÁLISE DE TENDÊNCIA (REGRESSÃO DE PRAIS-WINSTEN) ===

# A regressão de Prais-Winsten (PW) é usada para análise de tendência em
# séries temporais que possuem autocorrelação de primeira ordem (AR1).

# --- 8.1 Funções de Suporte (Durbin-Watson e Wrapper da Regressão) ---

#' Calcula a estatística de Durbin-Watson
#' @param resid Resíduos de um modelo linear
dw_stat <- function(resid){
  resid <- as.numeric(resid)
  if(length(resid) < 3 || any(!is.finite(resid))) return(NA_real_)
  sum(diff(resid)^2) / sum(resid^2)
}

#' Wrapper seguro para rodar a regressão de Prais-Winsten
#'
#' Tenta rodar `prais_winsten`. Se falhar (comum em séries curtas ou
#' com baixa variância), utiliza um fallback para `nlme::gls` com
#' correção AR(1), que é mais robusto.
#'
#' @param dat Dataframe para uma única região, deve conter `ddd_1000hab_dia` e `time_id`
#' @return Um tibble (1 linha) com as métricas do modelo (beta, p-valor, etc.)
run_pw_safe <- function(dat){
  dat <- dat %>% arrange(time_id)
  n <- nrow(dat)
  
  # Requer um N mínimo para a regressão
  if(n < 6) {
    return(tibble(
      beta=NA_real_, se=NA_real_, ci_low=NA_real_, ci_high=NA_real_,
      p_value=NA_real_, r2=NA_real_, n=n, rho=NA_real_,
      dw_original=NA_real_, dw_transformed=NA_real_, method="insufficient_n"
    ))
  }
  
  # 1. Roda OLS (linear simples) para obter o DW original
  fit_ols <- lm(ddd_1000hab_dia ~ time_id, data = dat)
  dw_orig <- dw_stat(residuals(fit_ols))
  
  # 2. Tenta rodar Prais-Winsten
  pw_try <- try(
    prais_winsten(ddd_1000hab_dia ~ time_id, data = dat, index = "time_id"),
    silent = TRUE
  )
  
  if(inherits(pw_try, "try-error")){
    # 3. FALLBACK: Se Prais falhar, tenta nlme::gls(corAR1)
    gls_try <- try({
      nlme::gls(ddd_1000hab_dia ~ time_id,
                data = dat,
                correlation = nlme::corAR1(form = ~ time_id))
    }, silent = TRUE)
    
    if(inherits(gls_try, "try-error")){
      # 4. FALHA TOTAL: Retorna NAs
      return(tibble(
        beta=NA_real_, se=NA_real_, ci_low=NA_real_, ci_high=NA_real_,
        p_value=NA_real_, r2=summary(fit_ols)$adj.r.squared, n=n,
        rho=NA_real_, dw_original=dw_orig, dw_transformed=NA_real_,
        method="failed_pw_and_gls"
      ))
    } else {
      # 5. SUCESSO (GLS): Extrai resultados do GLS
      sm <- summary(gls_try)
      co <- sm$tTable["time_id", ]
      se <- unname(co["Std.Error"])
      beta <- unname(co["Value"])
      df  <- n - 2
      t95 <- qt(0.975, df)
      return(tibble(
        beta=beta, se=se,
        ci_low=beta - t95*se, ci_high=beta + t95*se,
        p_value=unname(co["p-value"]),
        r2=summary(fit_ols)$adj.r.squared, # R² do OLS (GLS não tem R² direto)
        n=n,
        rho=as.numeric(nlme::coef(gls_try$modelStruct$corStruct, unconstrained = FALSE)),
        dw_original=dw_orig, dw_transformed=dw_stat(residuals(gls_try)),
        method="gls_ar1_fallback"
      ))
    }
  } else {
    # 6. SUCESSO (PRAIS): Extrai resultados do Prais-Winsten
    sm <- summary(pw_try)
    co <- sm$coefficients["time_id", ]
    se <- unname(co["Std. Error"])
    beta <- unname(co["Estimate"])
    df  <- n - 2
    t95 <- qt(0.975, df)
    return(tibble(
      beta=beta, se=se,
      ci_low=beta - t95*se, ci_high=beta + t95*se,
      p_value=unname(co["Pr(>|t|)"]),
      r2=sm$r.squared, n=n,
      rho=if(!is.null(pw_try$rho)) as.numeric(pw_try$rho) else NA_real_,
      dw_original=dw_orig, dw_transformed=dw_stat(residuals(pw_try)),
      method="prais_winsten"
    ))
  }
}

# --- 8.2 Preparação dos Dados para Análise ---
# Cria o índice temporal (time_id) sequencial para cada região
df_pw <- df_out %>%
  # Garante ordem e dados únicos
  arrange(regiao, ano, trimestre) %>%
  distinct(regiao, ano, trimestre, .keep_all = TRUE) %>%
  
  # Remove NAs/Infinitos que quebram a regressão
  filter(is.finite(ddd_1000hab_dia)) %>%
  
  group_by(regiao) %>%
  mutate(time_id = row_number()) %>% # Índice 1, 2, 3... por região
  ungroup()

# --- 8.3 Execução da Análise ---
# Roda a função 'run_pw_safe' para cada grupo (região)
res_pw <- df_pw %>%
  group_by(regiao) %>%
  group_modify(~ run_pw_safe(.x)) %>%
  ungroup()

# (Opcional) Conferir a coluna 'method' para ver qual modelo foi usado
# print(res_pw)


# === 9. FORMATAÇÃO E EXPORTAÇÃO DOS RESULTADOS ===

# --- 9.1 Função de Suporte (Formatar p-valor) ---
fmt_p <- function(p){
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

# --- 9.2 Criação da Tabela Final (Formato de Artigo) ---
tabela_artigo <- res_pw %>%
  mutate(
    # Calcula a inclinação (Beta) anual (Beta trimestral * 4)
    beta_year = beta * 4,
    ci_low_year = ci_low * 4,
    ci_high_year = ci_high * 4
  ) %>%
  # Formata as colunas para a tabela final
  transmute(
    Region = regiao,
    `Beta (per quarter)` = sprintf("%.4f", beta),
    `95% CI (per quarter)` = sprintf("[%.4f; %.4f]", ci_low, ci_high),
    `Beta (per year)` = sprintf("%.4f", beta_year),
    `95% CI (per year)` = sprintf("[%.4f; %.4f]", ci_low_year, ci_high_year),
    `p-value` = fmt_p(p_value),
    `Adj. R²` = sprintf("%.3f", r2),
    `AR(1) rho` = sprintf("%.3f", rho),
    `DW original` = sprintf("%.3f", dw_original),
    `DW transformed` = sprintf("%.3f", dw_transformed),
    method,
    n = n
  ) %>%
  arrange(Region)

# --- 9.3 Visualização e Exportação ---

# Garante que há apenas uma linha por região (caso haja duplicatas)
resultado_final <- tabela_artigo %>%
  group_by(Region) %>%
  summarise(across(everything(), ~ first(.)), .groups = "drop")

# Visualiza a tabela no console
print(resultado_final)

# Salva a tabela formatada em um arquivo .csv
write.csv(resultado_final, 
          "resultado_tendencia_risperidona.csv", 
          row.names = FALSE) # row.names = FALSE é importante!

# === FIM DO SCRIPT ===
