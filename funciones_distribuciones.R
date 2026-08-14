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

#jajajaj
mi_geometrico <- function(tipo, x, prob) {
  if (tipo == "d") {
    (1 - prob)^x * prob  # Lógica de masa de probabilidad (PMF)
  
    } else if (tipo == "p") {
      sum(
        sapply(0:x, function(t) {
          (1 - prob)^t * prob
        })
      ) # Lógica acumulada (CDF)
      
  } else if (tipo == "q") {
    limite <- max(100, ceiling(qgeom(0.9999, prob) + 10))
    pdf <- sapply(0:limite, function(t) {
      (1 - prob)^t * prob
    })
    which(cumsum(pdf) >= x)[1] - 1 # Lógica de cuantiles
    
  } else if (tipo == "r") {
    # Números aleatorios usando la transformada inversa o conteo de fallos
    replicate(x, {
      k <- 0
      while (runif(1) >= prob) {
        k <- k + 1
      }
      k
    })
  }
}

mi_hipergeometrico <- function(tipo, x, m, n, k) {
  # m: número de éxitos en la población, n: número de fracasos, k: número de extracciones
  if (tipo == "d") {
    # Función masa de probabilidad
    (choose(m, x) * choose(n, k - x)) / choose(m + n, k)
    
  } else if (tipo == "p") {
    # Función acumulada
    min_val <- max(0, k - n)
    max_val <- min(k, m)
    rango <- min_val:x
    sum(
      sapply(rango, function(t) {
        (choose(m, t) * choose(n, k - t)) / choose(m + n, k)
      })
    )
    
  } else if (tipo == "q") {
    # Cuantiles
    min_val <- max(0, k - n)
    max_val <- min(k, m)
    rango <- min_val:max_val
    pdf <- sapply(rango, function(t) {
      (choose(m, t) * choose(n, k - t)) / choose(m + n, k)
    })
    rango[which(cumsum(pdf) >= x)[1]]
    
  } else if (tipo == "r") {
    # Números aleatorios mediante simulación de urna con sample
    replicate(x, {
      urna <- rep(c(1, 0), c(m, n))
      extraccion <- sample(urna, k, replace = FALSE)
      sum(extraccion)
    })
  }
}

mi_F <- function(tipo, x, df1, df2) {
  if (tipo == "d") {
    # Función de densidad de la F de Snedecor
    num <- (df1 * df2)^(df1 / 2) * gamma((df1 + df2) / 2) * x^((df1 / 2) - 1)
    den <- gamma(df1 / 2) * gamma(df2 / 2) * (df2 + df1 * x)^((df1 + df2) / 2)
    num / den
    
  } else if (tipo == "p") {
    # Función acumulada
    integrate(function(t) {
      num <- (df1 * df2)^(df1 / 2) * gamma((df1 + df2) / 2) * t^((df1 / 2) - 1)
      den <- gamma(df1 / 2) * gamma(df2 / 2) * (df2 + df1 * t)^((df1 + df2) / 2)
      num / den
    }, 0, x)$value
    
  } else if (tipo == "q") {
    # Cuantiles
    limite_superior <- 1000
    uniroot(function(u) {
      integrate(function(t) {
        num <- (df1 * df2)^(df1 / 2) * gamma((df1 + df2) / 2) * t^((df1 / 2) - 1)
        den <- gamma(df1 / 2) * gamma(df2 / 2) * (df2 + df1 * t)^((df1 + df2) / 2)
        num / den
      }, 0, u)$value - x
    }, interval = c(0, limite_superior))$root
    
  } else if (tipo == "r") {
    # Una distribución F se forma como el cociente de dos Ji-cuadradas divididas entre sus grados de libertad
    replicate(x, {
      chi1 <- sum((sqrt(-2 * log(runif(df1))) * cos(2 * pi * runif(df1)))^2)
      chi2 <- sum((sqrt(-2 * log(runif(df2))) * cos(2 * pi * runif(df2)))^2)
      (chi1 / df1) / (chi2 / df2)
    })
  }
}

mi_gama <- function(tipo, x, shape, rate = 1) {
  # shape: parámetro de forma (alpha), rate: parámetro de tasa (beta)
  if (tipo == "d") {
    # Función de densidad
    (rate^shape / gamma(shape)) * x^(shape - 1) * exp(-rate * x)
    
  } else if (tipo == "p") {
    # Función acumulada
    integrate(function(t) {
      (rate^shape / gamma(shape)) * t^(shape - 1) * exp(-rate * t)
    }, 0, x)$value
    
  } else if (tipo == "q") {
    # Cuantiles
    limite_superior <- max(100, (shape / rate) + 10 * sqrt(shape / rate^2))
    uniroot(function(u) {
      integrate(function(t) {
        (rate^shape / gamma(shape)) * t^(shape - 1) * exp(-rate * t)
      }, 0, u)$value - x
    }, interval = c(0, limite_superior))$root
    
  } else if (tipo == "r") {
    # Método basado en la suma de exponenciales (para shape entero) o métodos generales
    # Usando el método de transformación/suma para shape entero o aproximación general:
    replicate(x, {
      # Si shape es entero, una Gamma(shape, rate) es la suma de 'shape' Exponenciales(rate)
      if (shape == floor(shape)) {
        sum(replicate(shape, -log(runif(1)) / rate))
      } else {
        # Aproximación general usando el método de aceptación-rechazo o sumas gamma
        sum(replicate(ceiling(shape), -log(runif(1)) / rate)) # Aproximación base
      }
    })
  }
}

mi_beta <- function(tipo, x, shape1, shape2) {
  # shape1: alpha, shape2: beta
  if (tipo == "d") {
    # Función de densidad
    beta_constante <- (gamma(shape1) * gamma(shape2)) / gamma(shape1 + shape2)
    (x^(shape1 - 1) * (1 - x)^(shape2 - 1)) / beta_constante
    
  } else if (tipo == "p") {
    # Función acumulada
    beta_constante <- (gamma(shape1) * gamma(shape2)) / gamma(shape1 + shape2)
    integrate(function(t) {
      (t^(shape1 - 1) * (1 - t)^(shape2 - 1)) / beta_constante
    }, 0, x)$value
    
  } else if (tipo == "q") {
    # Cuantiles
    uniroot(function(u) {
      beta_constante <- (gamma(shape1) * gamma(shape2)) / gamma(shape1 + shape2)
      integrate(function(t) {
        (t^(shape1 - 1) * (1 - t)^(shape2 - 1)) / beta_constante
      }, 0, u)$value - x
    }, interval = c(0, 1))$root
    
  } else if (tipo == "r") {
    # Una Beta(shape1, shape2) se obtiene como X / (X + Y) donde X ~ Gamma(shape1, 1) y Y ~ Gamma(shape2, 1)
    replicate(x, {
      g1 <- sum(-log(runif(ceiling(shape1))))
      g2 <- sum(-log(runif(ceiling(shape2))))
      g1 / (g1 + g2)
    })
  }
}


# (Mai agregará aquí el esqueleto de sus otras 4 distribuciones:
# mi_hipergeometrico, mi_F, mi_gama, mi_beta)