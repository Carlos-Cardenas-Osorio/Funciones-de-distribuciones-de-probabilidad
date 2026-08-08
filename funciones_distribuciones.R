# ==========================================
# PROYECTO: REPLICA DE DISTRIBUCIONES EN R
# ==========================================

# ------------------------------------------
# 1. FUNCIÓN MAESTRA (Despachador)
# ------------------------------------------
# Esta función recibe qué distribución y qué tipo (d, p, q, r) se desea.

simulador_dist <- function(distribucion, tipo, x, ...) {
  # distribucion: "binomial", "normal", "poisson", etc.
  # tipo: "d" (densidad), "p" (acumulada), "q" (cuantil), "r" (aleatorio)
  
  if (distribucion == "normal") {
    resultado <- normal_custom(tipo, x, ...)
  } else if (distribucion == "binomial") {
    resultado <- binomial_custom(tipo, x, ...)
  }
  # Agregar el resto con ifelse o switch...
  
  return(resultado)
}

# ------------------------------------------
# 2. SECCIÓN DE [TU NOMBRE]
# (Aquí van tus 6 distribuciones)
# ------------------------------------------

normal_custom <- function(tipo, x, mean = 0, sd = 1) {
  if (tipo == "d") {
    # Lógica de dnorm matemática aquí
  } else if (tipo == "p") {
    # Lógica de pnorm
  } # ...
}

# ------------------------------------------
# 3. SECCIÓN DE MAI
# (Aquí van las 5 distribuciones de Mai)
# ------------------------------------------

binomial_custom <- function(tipo, x, size, prob) {
  if (tipo == "d") {
    # Lógica de dbinom matemática aquí
  } # ...
}