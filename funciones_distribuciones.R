# =====================================================================
# PROYECTO: RÉPLICA DE DISTRIBUCIONES ESTADÍSTICAS EN R
# Objetivo: Replicar las funciones nativas (d, p, q, r) desde cero.
# =====================================================================

# ---------------------------------------------------------------------
# SECCIÓN 1: [Carlos] (6 Distribuciones)
# Distribuciones: Normal, Binomial, Poisson, Exponencial, t de Student, ji-cuadrado
# ---------------------------------------------------------------------

mi_normal <- function(tipo, x, mean = 0, sd = 1) {
  if (tipo == "d") {
    # Lógica de densidad (PDF)
  } else if (tipo == "p") {
    # Lógica acumulada (CDF)
  } else if (tipo == "q") {
    # Lógica de cuantiles
  } else if (tipo == "r") {
    # Lógica de números aleatorios (n = x)
  }
}

# (Agrega aquí el esqueleto de tus otras 5 distribuciones: 
# mi_binomial, mi_poisson, mi_exponencial, mi_tstudent, mi_jicuadrado)


# ---------------------------------------------------------------------
# SECCIÓN 2: MAI (5 Distribuciones)
# Distribuciones: Geométrico, Hipergeométrico, F, Gama, Beta
# ---------------------------------------------------------------------

mi_geometrico <- function(tipo, x, prob) {
  if (tipo == "d") {
    # Lógica de masa de probabilidad (PMF)
  } else if (tipo == "p") {
    # Lógica acumulada (CDF)
  } else if (tipo == "q") {
    # Lógica de cuantiles
  } else if (tipo == "r") {
    # Lógica de números aleatorios (n = x)
  }
}

# (Mai agregará aquí el esqueleto de sus otras 4 distribuciones:
# mi_hipergeometrico, mi_F, mi_gama, mi_beta)