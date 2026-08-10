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

mi_binomial <- function(tipo, x, size, prob) {
  if (tipo == "d") {
    # Función masa de probabilidad
    choose(size,x)*prob^x*(1-prob)^(size-x)
    
  } else if (tipo == "p") {
    # Función acumulada
    sum(
      sapply(0:x,function(t){choose(size,t)*prob^t*(1-prob)^(size-t)}
             )
      )
    
  } else if (tipo == "q") {
    #Acumulada 
    pdf <- sapply(0:size, 
            function(t){
           choose(size,t)*prob^t*(1-prob)^(size-t)
                  })
    #Calculo del cuantil
    which(cumsum(pdf)>=x)[1]-1
    
  } else if (tipo == "r") {
    # x numeros aleatorios
    replicate(x, 
              sum(runif(size)<prob)
              )
  }
}

mi_poisson <- function(tipo, x, lambda) {
  if (tipo == "d") {
    # Función masa de probabilidad
    (lambda^x * exp(-lambda)) / factorial(x)
    
  } else if (tipo == "p") {
    # Función acumulada
    sum(
      sapply(0:x, function(t) {
        (lambda^t * exp(-lambda)) / factorial(t)
      })
    )
    
  } else if (tipo == "q") {
    # Limite seguro
    limite <- max(100, ceiling(lambda + 10 * sqrt(lambda)))
    pdf <- sapply(0:limite, function(t) {
      (lambda^t * exp(-lambda)) / factorial(t)
    })
    which(cumsum(pdf) >= x)[1] - 1
    
  } else if (tipo == "r") {
    # Algoritmo clásico de multiplicación de uniformes
    replicate(x, {
      L <- exp(-lambda)
      k <- 0
      p <- 1
      while (p > L) {
        p <- p * runif(1)
        k <- k + 1
      }
      k - 1
    })
  }
}

mi_exponencial <- function(tipo, x, rate = 1) {
  if (tipo == "d") {
    # Función de densidad
    rate * exp(-rate * x)
    
  } else if (tipo == "p") {
    # Función acumulada 
    integrate(function(t) {
      rate * exp(-rate * t)
    }, 0, x)$value
    
  } else if (tipo == "q") {
    # Cuantiles
    uniroot(function(u) {
      integrate(function(t) {
        rate * exp(-rate * t)
      }, 0, u)$value - x
    }, interval = c(0, 1000 / rate))$root
    
  } else if (tipo == "r") {
    # Método de la transformada inversa
    replicate(x, -log(runif(1)) / rate)
  }
}

mi_tstudent <- function(tipo, x, df) {
  if (tipo == "d") {
    # Función de densidad
    numerador <- gamma((df + 1) / 2)
    denominador <- sqrt(df * pi) * gamma(df / 2)
    (numerador / denominador) * (1 + (x^2) / df)^(-(df + 1) / 2)
    
  } else if (tipo == "p") {
    # Función acumulada
    integrate(function(t) {
      numerador <- gamma((df + 1) / 2)
      denominador <- sqrt(df * pi) * gamma(df / 2)
      (numerador / denominador) * (1 + (t^2) / df)^(-(df + 1) / 2)
    }, -Inf, x)$value
    
  } else if (tipo == "q") {
    # Cuantiles
    uniroot(function(u) {
      integrate(function(t) {
        numerador <- gamma((df + 1) / 2)
        denominador <- sqrt(df * pi) * gamma(df / 2)
        (numerador / denominador) * (1 + (t^2) / df)^(-(df + 1) / 2)
      }, -Inf, u)$value - x
    }, interval = c(-1000, 1000))$root
    
  } else if (tipo == "r") {
    #T-Student se forma dividiendo una Normal Estándar 
    # entre la raíz de una Ji-cuadrado sobre sus grados de libertad
    replicate(x, {
      Z <- sqrt(-2 * log(runif(1))) * cos(2 * pi * runif(1)) # 1 Normal
      V <- sum((sqrt(-2 * log(runif(df))) * cos(2 * pi * runif(df)))^2) # 1 Ji-cuadrado
      Z / sqrt(V / df)
    })
  }
}

mi_jicuadrado <- function(tipo, x, df) {
  if (tipo == "d") {
    # Función de densidad
    (x^((df/2) - 1) * exp(-x/2)) / (2^(df/2) * gamma(df/2))
    
  } else if (tipo == "p") {
    # Función acumulada
    integrate(function(t) {
      (t^((df/2) - 1) * exp(-t/2)) / (2^(df/2) * gamma(df/2))
    }, 0, x)$value
    
  } else if (tipo == "q") {
    # Cuantiles
    limite_superior <- df + 100 * sqrt(2 * df) # Límite dinámico seguro
    uniroot(function(u) {
      integrate(function(t) {
        (t^((df/2) - 1) * exp(-t/2)) / (2^(df/2) * gamma(df/2))
      }, 0, u)$value - x
    }, interval = c(0, limite_superior))$root
    
  } else if (tipo == "r") {
    # Una Ji-cuadrado es la suma de 'df' Normales Estándar al cuadrado
    # Replicamos el método Box-Muller de tu función Normal y lo elevamos al cuadrado
    replicate(x, {
      normales_estandar <- sqrt(-2 * log(runif(df))) * cos(2 * pi * runif(df))
      sum(normales_estandar^2)
    })
  }
}


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