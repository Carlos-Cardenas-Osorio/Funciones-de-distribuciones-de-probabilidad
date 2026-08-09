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
    #Función de densidad
    1/(sd*sqrt(2*pi))*exp(-1/2*((x-mean)/sd)^2)
    
  } else if (tipo == "p") {
    #Función acumulada
    integrate( function(t){ 
      1/(sd*sqrt(2*pi))*exp(-1/2*((t-mean)/sd)^2)
      },-Inf,x)$value
    
  } else if (tipo == "q") {
    #Raíz
    uniroot(function(u){
      #Función acumulada
      integrate( function(t){    
        1/(sd*sqrt(2*pi))*exp(-1/2*((t-mean)/sd)^2)
      },-Inf,u)$value -x
    },interval = c(mean-10*sd,mean+10*sd))$root
    
  } else if (tipo == "r") {
    # Transformada de Box-Muller
    U_1 <- runif(x)
    U_2 <- runif(x)
    Z<-sqrt(-2*log(U_1))*cos(2*pi*U_2)
    X<-Z*sd+mean
    return(X)
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