# =====================================================================
# SISTEMA INTERACTIVO DE PROBABILIDAD Y DECISIONES
# =====================================================================

menu_principal <- function() {
  repeat {
    cat("\n========================================\n")
    cat("   MENÚ PRINCIPAL - HERRAMIENTAS\n")
    cat("========================================\n")
    cat("1. Calcular Distribuciones de Probabilidad\n")
    cat("2. Teoría de Decisiones bajo Incertidumbre\n")
    cat("3. Salir\n")
    cat("========================================\n")
    
    opcion <- readline(prompt = "Elige una opción (1-3): ")
    
    switch(opcion,
           "1" = {
             cat("\n[Iniciando módulo de Distribuciones...]\n")
             menu_distribuciones()
           },
           "2" = {
             cat("\n[Iniciando módulo de Teoría de Decisiones...]\n")
             menu_decisiones()
           },
           "3" = {
             cat("\n¡Hasta luego! Cerrando el sistema...\n")
             break # Rompe el ciclo repeat y termina la ejecución
           },
           {
             # Opción por defecto si el usuario teclea algo incorrecto
             cat("\n❌ Opción no válida. Por favor, ingresa 1, 2 o 3.\n")
           }
    )
  }
}


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
    # Usando logaritmos para evitar Inf / Inf
    exp(x * log(lambda) - lambda - lfactorial(x))
    
  } else if (tipo == "p") {
    # Función acumulada con logaritmos
    sum(
      sapply(0:x, function(t) {
        exp(t * log(lambda) - lambda - lfactorial(t))
      })
    )
    
  } else if (tipo == "q") {
    # Límite seguro con logaritmos
    limite <- max(100, ceiling(lambda + 10 * sqrt(lambda)))
    pdf <- sapply(0:limite, function(t) {
      exp(t * log(lambda) - lambda - lfactorial(t))
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

# =====================================================================
# FASE 2: MÓDULO DE DISTRIBUCIONES (COMPLETO)
# =====================================================================

menu_distribuciones <- function() {
  cat("\n--- SELECCIÓN DE DISTRIBUCIÓN ---\n")
  cat("1. Normal\n2. Binomial\n3. Poisson\n4. Exponencial\n")
  cat("5. t de Student\n6. Ji-cuadrado\n7. Geométrica\n8. Hipergeométrica\n")
  cat("9. F de Snedecor\n10. Gamma\n11. Beta\n")
  cat("0. Volver al menú principal\n")
  
  dist_opcion <- readline(prompt = "Elige la distribución (0-11): ")
  
  if (dist_opcion == "0") return()
  
  params <- list()
  nombre_dist <- ""
  
  # 1. PEDIR PARÁMETROS SEGÚN LA DISTRIBUCIÓN
  if (dist_opcion == "1") {
    nombre_dist <- "Normal"
    params$mean <- as.numeric(readline(prompt = "Ingresa la media (mean): "))
    params$sd <- as.numeric(readline(prompt = "Ingresa la desviación estándar (sd): "))
  } else if (dist_opcion == "2") {
    nombre_dist <- "Binomial"
    params$size <- as.numeric(readline(prompt = "Ingresa el número de ensayos (size): "))
    params$prob <- as.numeric(readline(prompt = "Ingresa la probabilidad de éxito (prob): "))
  } else if (dist_opcion == "3") {
    nombre_dist <- "Poisson"
    params$lambda <- as.numeric(readline(prompt = "Ingresa la tasa promedio (lambda): "))
  } else if (dist_opcion == "4") {
    nombre_dist <- "Exponencial"
    params$rate <- as.numeric(readline(prompt = "Ingresa la tasa (rate = 1/lambda): "))
  } else if (dist_opcion == "5") {
    nombre_dist <- "t de Student"
    params$df <- as.numeric(readline(prompt = "Ingresa los grados de libertad (df): "))
  } else if (dist_opcion == "6") {
    nombre_dist <- "Ji-cuadrado"
    params$df <- as.numeric(readline(prompt = "Ingresa los grados de libertad (df): "))
  } else if (dist_opcion == "7") {
    nombre_dist <- "Geométrica"
    params$prob <- as.numeric(readline(prompt = "Ingresa la probabilidad de éxito (prob): "))
  } else if (dist_opcion == "8") {
    nombre_dist <- "Hipergeométrica"
    params$m <- as.numeric(readline(prompt = "Éxitos en la población (m): "))
    params$n <- as.numeric(readline(prompt = "Fracasos en la población (n): "))
    params$k <- as.numeric(readline(prompt = "Número de extracciones (k): "))
  } else if (dist_opcion == "9") {
    nombre_dist <- "F de Snedecor"
    params$df1 <- as.numeric(readline(prompt = "Grados de libertad del numerador (df1): "))
    params$df2 <- as.numeric(readline(prompt = "Grados de libertad del denominador (df2): "))
  } else if (dist_opcion == "10") {
    nombre_dist <- "Gamma"
    params$shape <- as.numeric(readline(prompt = "Parámetro de forma (shape / alpha): "))
    params$rate <- as.numeric(readline(prompt = "Parámetro de tasa (rate / beta): "))
  } else if (dist_opcion == "11") {
    nombre_dist <- "Beta"
    params$shape1 <- as.numeric(readline(prompt = "Parámetro shape1 (alpha): "))
    params$shape2 <- as.numeric(readline(prompt = "Parámetro shape2 (beta): "))
  } else {
    cat("\n❌ Opción no válida.\n")
    return()
  }
  
  # 2. PEDIR QUÉ SE DESEA CALCULAR
  cat(paste0("\n--- DISTRIBUCIÓN ", toupper(nombre_dist), " ---\n"))
  cat("1. Función de Densidad / Masa f(x)\n")
  cat("2. Función Acumulada F(x)\n")
  cat("3. Función de Supervivencia S(x) = 1 - F(x)\n")
  cat("4. Cuantil (Dado un porcentaje)\n")
  cat("5. Generar números aleatorios\n")
  cat("6. Mostrar Media, Varianza y Desviación Estándar teóricas\n")
  
  calc_opcion <- readline(prompt = "¿Qué deseas calcular? (1-6): ")
  
  # 3. PEDIR EL VALOR DE EVALUACIÓN 'X' (Si no es la opción 6)
  x_val <- NA
  if (calc_opcion %in% c("1", "2", "3", "4", "5")) {
    if (calc_opcion == "4") {
      x_val <- as.numeric(readline(prompt = "Ingresa la probabilidad objetivo (ej. 0.95): "))
    } else if (calc_opcion == "5") {
      x_val <- as.numeric(readline(prompt = "Ingresa la cantidad de números a generar (n): "))
    } else {
      x_val <- as.numeric(readline(prompt = "Ingresa el valor de evaluación (x): "))
    }
  }
  
  # 4. EJECUTAR EL CÁLCULO SELECCIONADO
  cat("\n--- RESULTADO ---\n")
  
  # Bloques de cálculo para las 11 distribuciones
  if (dist_opcion == "1") { # NORMAL
    if (calc_opcion == "1") cat("f(x) =", mi_normal("d", x_val, params$mean, params$sd), "\n")
    if (calc_opcion == "2") cat("F(x) =", mi_normal("p", x_val, params$mean, params$sd), "\n")
    if (calc_opcion == "3") cat("S(x) =", 1 - mi_normal("p", x_val, params$mean, params$sd), "\n")
    if (calc_opcion == "4") cat("Cuantil =", mi_normal("q", x_val, params$mean, params$sd), "\n")
    if (calc_opcion == "5") { cat("Aleatorios:\n"); print(mi_normal("r", x_val, params$mean, params$sd)) }
    if (calc_opcion == "6") {
      cat("Media E[X] =", params$mean, "\n")
      cat("Varianza Var(X) =", params$sd^2, "\n")
      cat("Desviación Estándar =", params$sd, "\n")
    }
    
  } else if (dist_opcion == "2") { # BINOMIAL
    if (calc_opcion == "1") cat("f(x) =", mi_binomial("d", x_val, params$size, params$prob), "\n")
    if (calc_opcion == "2") cat("F(x) =", mi_binomial("p", x_val, params$size, params$prob), "\n")
    if (calc_opcion == "3") cat("S(x) =", 1 - mi_binomial("p", x_val, params$size, params$prob), "\n")
    if (calc_opcion == "4") cat("Cuantil =", mi_binomial("q", x_val, params$size, params$prob), "\n")
    if (calc_opcion == "5") { cat("Aleatorios:\n"); print(mi_binomial("r", x_val, params$size, params$prob)) }
    if (calc_opcion == "6") {
      cat("Media E[X] =", params$size * params$prob, "\n")
      cat("Varianza Var(X) =", params$size * params$prob * (1 - params$prob), "\n")
      cat("Desviación Estándar =", sqrt(params$size * params$prob * (1 - params$prob)), "\n")
    }
    
  } else if (dist_opcion == "3") { # POISSON
    if (calc_opcion == "1") cat("f(x) =", mi_poisson("d", x_val, params$lambda), "\n")
    if (calc_opcion == "2") cat("F(x) =", mi_poisson("p", x_val, params$lambda), "\n")
    if (calc_opcion == "3") cat("S(x) =", 1 - mi_poisson("p", x_val, params$lambda), "\n")
    if (calc_opcion == "4") cat("Cuantil =", mi_poisson("q", x_val, params$lambda), "\n")
    if (calc_opcion == "5") { cat("Aleatorios:\n"); print(mi_poisson("r", x_val, params$lambda)) }
    if (calc_opcion == "6") {
      cat("Media E[X] =", params$lambda, "\n")
      cat("Varianza Var(X) =", params$lambda, "\n")
      cat("Desviación Estándar =", sqrt(params$lambda), "\n")
    }
    
  } else if (dist_opcion == "4") { # EXPONENCIAL
    if (calc_opcion == "1") cat("f(x) =", mi_exponencial("d", x_val, params$rate), "\n")
    if (calc_opcion == "2") cat("F(x) =", mi_exponencial("p", x_val, params$rate), "\n")
    if (calc_opcion == "3") cat("S(x) =", 1 - mi_exponencial("p", x_val, params$rate), "\n")
    if (calc_opcion == "4") cat("Cuantil =", mi_exponencial("q", x_val, params$rate), "\n")
    if (calc_opcion == "5") { cat("Aleatorios:\n"); print(mi_exponencial("r", x_val, params$rate)) }
    if (calc_opcion == "6") {
      cat("Media E[X] =", 1 / params$rate, "\n")
      cat("Varianza Var(X) =", 1 / (params$rate^2), "\n")
      cat("Desviación Estándar =", sqrt(1 / (params$rate^2)), "\n")
    }
    
  } else if (dist_opcion == "5") { # T-STUDENT
    if (calc_opcion == "1") cat("f(x) =", mi_tstudent("d", x_val, params$df), "\n")
    if (calc_opcion == "2") cat("F(x) =", mi_tstudent("p", x_val, params$df), "\n")
    if (calc_opcion == "3") cat("S(x) =", 1 - mi_tstudent("p", x_val, params$df), "\n")
    if (calc_opcion == "4") cat("Cuantil =", mi_tstudent("q", x_val, params$df), "\n")
    if (calc_opcion == "5") { cat("Aleatorios:\n"); print(mi_tstudent("r", x_val, params$df)) }
    if (calc_opcion == "6") {
      cat("Media E[X] =", ifelse(params$df > 1, 0, "Indefinida (df <= 1)"), "\n")
      var_t <- ifelse(params$df > 2, params$df / (params$df - 2), ifelse(params$df > 1, Inf, "Indefinida"))
      cat("Varianza Var(X) =", var_t, "\n")
      cat("Desviación Estándar =", if(is.numeric(var_t)) sqrt(var_t) else "Indefinida", "\n")
    }
    
  } else if (dist_opcion == "6") { # JI-CUADRADO
    if (calc_opcion == "1") cat("f(x) =", mi_jicuadrado("d", x_val, params$df), "\n")
    if (calc_opcion == "2") cat("F(x) =", mi_jicuadrado("p", x_val, params$df), "\n")
    if (calc_opcion == "3") cat("S(x) =", 1 - mi_jicuadrado("p", x_val, params$df), "\n")
    if (calc_opcion == "4") cat("Cuantil =", mi_jicuadrado("q", x_val, params$df), "\n")
    if (calc_opcion == "5") { cat("Aleatorios:\n"); print(mi_jicuadrado("r", x_val, params$df)) }
    if (calc_opcion == "6") {
      cat("Media E[X] =", params$df, "\n")
      cat("Varianza Var(X) =", 2 * params$df, "\n")
      cat("Desviación Estándar =", sqrt(2 * params$df), "\n")
    }
    
  } else if (dist_opcion == "7") { # GEOMÉTRICA
    if (calc_opcion == "1") cat("f(x) =", mi_geometrico("d", x_val, params$prob), "\n")
    if (calc_opcion == "2") cat("F(x) =", mi_geometrico("p", x_val, params$prob), "\n")
    if (calc_opcion == "3") cat("S(x) =", 1 - mi_geometrico("p", x_val, params$prob), "\n")
    if (calc_opcion == "4") cat("Cuantil =", mi_geometrico("q", x_val, params$prob), "\n")
    if (calc_opcion == "5") { cat("Aleatorios:\n"); print(mi_geometrico("r", x_val, params$prob)) }
    if (calc_opcion == "6") {
      cat("Media E[X] =", (1 - params$prob) / params$prob, "\n")
      cat("Varianza Var(X) =", (1 - params$prob) / (params$prob^2), "\n")
      cat("Desviación Estándar =", sqrt((1 - params$prob) / (params$prob^2)), "\n")
    }
    
  } else if (dist_opcion == "8") { # HIPERGEOMÉTRICA
    if (calc_opcion == "1") cat("f(x) =", mi_hipergeometrico("d", x_val, params$m, params$n, params$k), "\n")
    if (calc_opcion == "2") cat("F(x) =", mi_hipergeometrico("p", x_val, params$m, params$n, params$k), "\n")
    if (calc_opcion == "3") cat("S(x) =", 1 - mi_hipergeometrico("p", x_val, params$m, params$n, params$k), "\n")
    if (calc_opcion == "4") cat("Cuantil =", mi_hipergeometrico("q", x_val, params$m, params$n, params$k), "\n")
    if (calc_opcion == "5") { cat("Aleatorios:\n"); print(mi_hipergeometrico("r", x_val, params$m, params$n, params$k)) }
    if (calc_opcion == "6") {
      N_pob <- params$m + params$n
      p_exito <- params$m / N_pob
      media_hip <- params$k * p_exito
      var_hip <- params$k * p_exito * (1 - p_exito) * ((N_pob - params$k) / (N_pob - 1))
      cat("Media E[X] =", media_hip, "\n")
      cat("Varianza Var(X) =", var_hip, "\n")
      cat("Desviación Estándar =", sqrt(var_hip), "\n")
    }
    
  } else if (dist_opcion == "9") { # F DE SNEDECOR
    if (calc_opcion == "1") cat("f(x) =", mi_F("d", x_val, params$df1, params$df2), "\n")
    if (calc_opcion == "2") cat("F(x) =", mi_F("p", x_val, params$df1, params$df2), "\n")
    if (calc_opcion == "3") cat("S(x) =", 1 - mi_F("p", x_val, params$df1, params$df2), "\n")
    if (calc_opcion == "4") cat("Cuantil =", mi_F("q", x_val, params$df1, params$df2), "\n")
    if (calc_opcion == "5") { cat("Aleatorios:\n"); print(mi_F("r", x_val, params$df1, params$df2)) }
    if (calc_opcion == "6") {
      cat("Media E[X] =", ifelse(params$df2 > 2, params$df2 / (params$df2 - 2), "Indefinida (df2 <= 2)"), "\n")
      if (params$df2 > 4) {
        num_f <- 2 * (params$df2^2) * (params$df1 + params$df2 - 2)
        den_f <- params$df1 * ((params$df2 - 2)^2) * (params$df2 - 4)
        var_f <- num_f / den_f
        cat("Varianza Var(X) =", var_f, "\n")
        cat("Desviación Estándar =", sqrt(var_f), "\n")
      } else {
        cat("Varianza Var(X) = Indefinida (df2 <= 4)\n")
        cat("Desviación Estándar = Indefinida\n")
      }
    }
    
  } else if (dist_opcion == "10") { # GAMMA
    if (calc_opcion == "1") cat("f(x) =", mi_gama("d", x_val, params$shape, params$rate), "\n")
    if (calc_opcion == "2") cat("F(x) =", mi_gama("p", x_val, params$shape, params$rate), "\n")
    if (calc_opcion == "3") cat("S(x) =", 1 - mi_gama("p", x_val, params$shape, params$rate), "\n")
    if (calc_opcion == "4") cat("Cuantil =", mi_gama("q", x_val, params$shape, params$rate), "\n")
    if (calc_opcion == "5") { cat("Aleatorios:\n"); print(mi_gama("r", x_val, params$shape, params$rate)) }
    if (calc_opcion == "6") {
      cat("Media E[X] =", params$shape / params$rate, "\n")
      cat("Varianza Var(X) =", params$shape / (params$rate^2), "\n")
      cat("Desviación Estándar =", sqrt(params$shape / (params$rate^2)), "\n")
    }
    
  } else if (dist_opcion == "11") { # BETA
    if (calc_opcion == "1") cat("f(x) =", mi_beta("d", x_val, params$shape1, params$shape2), "\n")
    if (calc_opcion == "2") cat("F(x) =", mi_beta("p", x_val, params$shape1, params$shape2), "\n")
    if (calc_opcion == "3") cat("S(x) =", 1 - mi_beta("p", x_val, params$shape1, params$shape2), "\n")
    if (calc_opcion == "4") cat("Cuantil =", mi_beta("q", x_val, params$shape1, params$shape2), "\n")
    if (calc_opcion == "5") { cat("Aleatorios:\n"); print(mi_beta("r", x_val, params$shape1, params$shape2)) }
    if (calc_opcion == "6") {
      suma_ab <- params$shape1 + params$shape2
      media_beta <- params$shape1 / suma_ab
      var_beta <- (params$shape1 * params$shape2) / ((suma_ab^2) * (suma_ab + 1))
      cat("Media E[X] =", media_beta, "\n")
      cat("Varianza Var(X) =", var_beta, "\n")
      cat("Desviación Estándar =", sqrt(var_beta), "\n")
    }
  }
}

# =====================================================================
# FASE 3: MÓDULO DE TEORÍA DE DECISIONES (SIN HURWICZ)
# =====================================================================

menu_decisiones <- function() {
  cat("\n--- MÓDULO DE TEORÍA DE DECISIONES ---\n")
  
  n_alt <- as.integer(readline(prompt = "¿Cuántas alternativas (filas) tienes? "))
  n_est <- as.integer(readline(prompt = "¿Cuántos estados de la naturaleza (columnas) tienes? "))
  
  if (is.na(n_alt) || is.na(n_est) || n_alt < 1 || n_est < 1) {
    cat("\n❌ Entradas no válidas. Deben ser números enteros mayores a 0.\n")
    return()
  }
  
  # Inicializamos la matriz vacía
  matriz <- matrix(0, nrow = n_alt, ncol = n_est)
  rownames(matriz) <- paste0("A", 1:n_alt)
  colnames(matriz) <- paste0("N", 1:n_est)
  
  cat("\n[Ingresando valores de la matriz]\n")
  for (i in 1:n_alt) {
    for (j in 1:n_est) {
      val <- as.numeric(readline(prompt = sprintf("Valor para Alternativa %d (A%d), Estado %d (N%d): ", i, i, j, j)))
      matriz[i, j] <- val
    }
  }
  
  cat("\nLa matriz ingresada es:\n")
  print(matriz)
  
  cat("\n¿Qué tipo de problema es?\n")
  cat("1. Beneficios / Ganancias (Maximizar)\n")
  cat("2. Costos / Tiempos (Minimizar)\n")
  tipo_prob <- readline(prompt = "Elige una opción (1-2): ")
  
  if (!(tipo_prob %in% c("1", "2"))) {
    cat("\n❌ Opción no válida.\n")
    return()
  }
  
  cat("\n¿Deseas incluir probabilidades para calcular el Valor Monetario Esperado (VME)?\n")
  cat("1. Sí (Riesgo)\n")
  cat("2. No, solo criterios sin probabilidades (Incertidumbre)\n")
  usa_probs <- readline(prompt = "Elige una opción (1-2): ")
  
  probs <- NULL
  if (usa_probs == "1") {
    probs <- numeric(n_est)
    cat("\n[Ingresando probabilidades para los estados de la naturaleza]\n")
    for (j in 1:n_est) {
      probs[j] <- as.numeric(readline(prompt = sprintf("Probabilidad para el Estado %d (N%d): ", j, j)))
    }
    # Pequeña validación por si la suma no es 1
    if (abs(sum(probs) - 1) > 0.001) {
      cat("\n⚠️ Advertencia: Las probabilidades no suman 1. Suma actual:", sum(probs), "\n")
    }
  }
  
  cat("\n========================================\n")
  cat("          RESULTADOS DEL ANÁLISIS       \n")
  cat("========================================\n")
  
  nombres_alt <- rownames(matriz)
  
  if (tipo_prob == "1") { # LÓGICA DE BENEFICIOS
    opt <- apply(matriz, 1, max)
    pes <- apply(matriz, 1, min)
    lap <- apply(matriz, 1, mean)
    
    # Matriz de arrepentimiento de Savage (Beneficios)
    max_cols <- apply(matriz, 2, max)
    matriz_arr <- sweep(matriz, 2, max_cols, "-") * -1
    sav <- apply(matriz_arr, 1, max)
    
    mejor_opt <- nombres_alt[which.max(opt)]
    mejor_pes <- nombres_alt[which.max(pes)]
    mejor_lap <- nombres_alt[which.max(lap)]
    mejor_sav <- nombres_alt[which.min(sav)] # El arrepentimiento siempre se minimiza
    
  } else { # LÓGICA DE COSTOS
    opt <- apply(matriz, 1, min)
    pes <- apply(matriz, 1, max)
    lap <- apply(matriz, 1, mean)
    
    # Matriz de arrepentimiento de Savage (Costos)
    min_cols <- apply(matriz, 2, min)
    matriz_arr <- sweep(matriz, 2, min_cols, "-")
    sav <- apply(matriz_arr, 1, max)
    
    mejor_opt <- nombres_alt[which.min(opt)]
    mejor_pes <- nombres_alt[which.min(pes)]
    mejor_lap <- nombres_alt[which.min(lap)]
    mejor_sav <- nombres_alt[which.min(sav)]
  }
  
  # Agrupando resultados
  resultados <- data.frame(Optimista = opt, Pesimista = pes, Laplace = lap, Max_Arrepentimiento = sav)
  
  if (usa_probs == "1") {
    vme <- as.vector(matriz %*% probs)
    resultados$VME <- vme
    mejor_vme <- if(tipo_prob == "1") nombres_alt[which.max(vme)] else nombres_alt[which.min(vme)]
  }
  
  print(resultados)
  
  cat("\n--- MATRIZ DE ARREPENTIMIENTO (SAVAGE) ---\n")
  print(matriz_arr)
  
  cat("\n--- MEJORES ALTERNATIVAS POR CRITERIO ---\n")
  cat("Optimista:               ", mejor_opt, "\n")
  cat("Pesimista:               ", mejor_pes, "\n")
  cat("Laplace (Promedio):      ", mejor_lap, "\n")
  cat("Savage (Arrepentimiento):", mejor_sav, "\n")
  if (usa_probs == "1") {
    cat("VME (Valor Mon. Esp.):   ", mejor_vme, "\n")
  }
  cat("========================================\n")
}
