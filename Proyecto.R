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
    replicate(x, {
      g1 <- sum(-log(runif(ceiling(shape1))))
      g2 <- sum(-log(runif(ceiling(shape2))))
      g1 / (g1 + g2)
    })
  }
}

# =====================================================================
# FASE 2: MÓDULO DE DISTRIBUCIONES (REESTRUCTURADO LÓGICAMENTE)
# =====================================================================

menu_distribuciones <- function() {
  cat("\n--- SELECCIÓN DE TIPO DE VARIABLE ---\n")
  cat("1. Distribuciones Discretas\n")
  cat("2. Distribuciones Continuas\n")
  cat("0. Volver al menú principal\n")
  
  tipo_var <- readline(prompt = "Elige una opción (0-2): ")
  if (tipo_var == "0") return()
  
  cat("\n--- TIPO DE INFORMACIÓN ---\n")
  cat("1. Conozco los parámetros teóricos\n")
  cat("2. Desconozco los parámetros (Ingresar muestra empírica para comparar modelos)\n")
  
  conocimiento <- readline(prompt = "Elige una opción (1-2): ")
  
  if (conocimiento == "1") {
    # RUTA A: PARÁMETROS CONOCIDOS
    if (tipo_var == "1") {
      cat("\n--- DISCRETAS ---\n1. Binomial\n2. Poisson\n3. Geométrica\n4. Hipergeométrica\n")
      opc <- readline(prompt = "Elige (1-4): ")
      dist_elegida <- switch(opc, "1"="binomial", "2"="poisson", "3"="geometrica", "4"="hipergeometrica", "error")
    } else {
      cat("\n--- CONTINUAS ---\n1. Normal\n2. Exponencial\n3. t de Student\n4. Ji-cuadrado\n5. F\n6. Gamma\n7. Beta\n")
      opc <- readline(prompt = "Elige (1-7): ")
      dist_elegida <- switch(opc, "1"="normal", "2"="exponencial", "3"="tstudent", "4"="jicuadrado", "5"="f", "6"="gamma", "7"="beta", "error")
    }
    
    if (dist_elegida == "error") return(cat("\n❌ Opción no válida.\n"))
    
    params <- capturar_parametros_manual(dist_elegida)
    menu_probabilidades_avanzado(dist_elegida, params, tipo_var)
    
  } else if (conocimiento == "2") {
    # RUTA B: DATOS EMPÍRICOS (PANEL COMPARATIVO)
    menu_ajuste_comparativo(tipo_var)
    
  } else {
    cat("\n❌ Opción no válida.\n")
  }
}

# =====================================================================
# FASE 3: MÓDULO DE TEORÍA DE DECISIONES
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


# =====================================================================
# PASO 2: ESTIMACIÓN (MLE), CORRELACIÓN Y AJUSTE VISUAL
# =====================================================================

ajuste_maxima_verosimilitud <- function(dist_elegida, tipo_var) {
  cat("\n--- INGRESO DE DATOS (MUESTRA) ---\n")
  cat("Ingresa los datos de tu muestra de X separados por comas (ej. 1.2, 3.4, 2.1):\n")
  
  datos_x_str <- readline(prompt = "Datos X: ")
  x <- as.numeric(unlist(strsplit(datos_x_str, ",")))
  x <- x[!is.na(x)] # Limpieza por si hay espacios extra
  
  if (length(x) == 0) {
    cat("\n❌ Error: No se detectaron datos válidos.\n")
    return(NULL)
  }
  
  cat(sprintf("\n✅ Se cargaron %d datos exitosamente.\n", length(x)))
  
  # ---------------------------------------------------------
  # COVARIANZA Y CORRELACIÓN
  # ---------------------------------------------------------
  cat("\n¿Deseas ingresar una segunda variable (Y) para calcular covarianza y correlación?\n")
  cat("1. Sí\n2. No\n")
  opc_cor <- readline(prompt = "Elige (1-2): ")
  
  if (opc_cor == "1") {
    cat("Ingresa los datos de Y separados por comas:\n")
    datos_y_str <- readline(prompt = "Datos Y: ")
    y <- as.numeric(unlist(strsplit(datos_y_str, ",")))
    y <- y[!is.na(y)]
    
    if (length(x) == length(y)) {
      cat("\n--- ESTADÍSTICA BIVARIADA ---\n")
      cat("Covarianza Cov(X,Y):", cov(x, y), "\n")
      cat("Correlación Cor(X,Y):", cor(x, y), "\n")
      cat("-----------------------------\n")
    } else {
      cat("\n Advertencia: Los vectores X e Y no tienen el mismo tamaño. Se omitirá el cálculo de correlación.\n")
    }
  }
  
  # ---------------------------------------------------------
  # ESTIMACIÓN DE PARÁMETROS (MLE / Momentos)
  # ---------------------------------------------------------
  params <- list()
  cat("\n[Calculando estimadores para el ajuste...]\n")
  
  if (dist_elegida == "normal") {
    params$mean <- mean(x)
    params$sd <- sd(x) 
    cat(sprintf("Parámetros estimados -> Media: %.4f, Desv. Est: %.4f\n", params$mean, params$sd))
    
  } else if (dist_elegida == "exponencial") {
    params$rate <- 1 / mean(x)
    cat(sprintf("Parámetro estimado -> Rate (lambda): %.4f\n", params$rate))
    
  } else if (dist_elegida == "poisson") {
    params$lambda <- mean(x)
    cat(sprintf("Parámetro estimado -> Lambda: %.4f\n", params$lambda))
    
  } else if (dist_elegida == "binomial") {
    params$size <- as.numeric(readline(prompt = "Para Binomial, necesitamos el número total de ensayos por observación (size): "))
    params$prob <- mean(x) / params$size
    cat(sprintf("Parámetro estimado -> Probabilidad (p): %.4f\n", params$prob))
    
  } else if (dist_elegida == "geometrica") {
    params$prob <- 1 / (1 + mean(x))
    cat(sprintf("Parámetro estimado -> Probabilidad (p): %.4f\n", params$prob))
    
  } else if (dist_elegida == "gamma") {
    # Estimación por Método de Momentos para Gamma
    varianza_muestral <- var(x)
    params$shape <- (mean(x)^2) / varianza_muestral
    params$rate <- mean(x) / varianza_muestral
    cat(sprintf("Parámetros estimados -> Shape: %.4f, Rate: %.4f\n", params$shape, params$rate))
    
  } else {
    cat("\n Para esta distribución compleja se recomienda usar métodos numéricos (ej. MASS::fitdistr) o parámetros conocidos.\n")
    return(NULL)
  }
  
  # ---------------------------------------------------------
  # AJUSTE VISUAL (GRÁFICA)
  # ---------------------------------------------------------
  cat("\nGenerando gráfica de ajuste visual en la ventana de RStudio (Plots)...\n")
  
  if (tipo_var == "2") { # Ajuste para variables continuas
    hist(x, probability = TRUE, col = "lightblue", border = "white",
         main = paste("Ajuste Teórico vs Datos:", toupper(dist_elegida)), 
         xlab = "Valores de X", ylab = "Densidad")
    
    # Secuencia de valores para dibujar la curva suavemente
    x_seq <- seq(min(x), max(x), length.out = 100)
    
    # Evaluar la función teórica "hecha a mano"
    if (dist_elegida == "normal") {
      y_teorico <- sapply(x_seq, function(v) mi_normal("d", v, params$mean, params$sd))
    } else if (dist_elegida == "exponencial") {
      y_teorico <- sapply(x_seq, function(v) mi_exponencial("d", v, params$rate))
    } else if (dist_elegida == "gamma") {
      y_teorico <- sapply(x_seq, function(v) mi_gama("d", v, params$shape, params$rate))
    }
    
    lines(x_seq, y_teorico, col = "red", lwd = 2)
    legend("topright", legend=c("Datos Muestrales", "Modelo Teórico"), fill=c("lightblue", NA), 
           border=c("black", NA), col=c(NA, "red"), lwd=c(NA, 2), bty = "n")
    
  } else { # Ajuste para variables discretas
    freq_relativa <- table(x) / length(x)
    plot(freq_relativa, type = "h", lwd = 10, col = "lightblue",
         main = paste("Ajuste Teórico vs Datos:", toupper(dist_elegida)),
         xlab = "Valores de X", ylab = "Masa de Probabilidad")
    
    x_unicos <- as.numeric(names(freq_relativa))
    
    if (dist_elegida == "poisson") {
      y_teorico <- sapply(x_unicos, function(v) mi_poisson("d", v, params$lambda))
    } else if (dist_elegida == "binomial") {
      y_teorico <- sapply(x_unicos, function(v) mi_binomial("d", v, params$size, params$prob))
    } else if (dist_elegida == "geometrica") {
      y_teorico <- sapply(x_unicos, function(v) mi_geometrico("d", v, params$prob))
    }
    
    points(x_unicos, y_teorico, col = "red", pch = 16, cex = 1.5)
    legend("topright", legend=c("Datos Muestrales", "Modelo Teórico"), pch=c(15, 16), 
           col=c("lightblue", "red"), bty = "n")
  }
  
  return(params) # Retornamos los parámetros para usarlos en el menú de probabilidades
}

# =====================================================================
# PASO 3 Y 4: MENÚ AVANZADO DE PROBABILIDADES Y CAPTURA MANUAL
# =====================================================================

capturar_parametros_manual <- function(dist_elegida) {
  params <- list()
  if (dist_elegida == "normal") {
    params$mean <- as.numeric(readline(prompt = "Ingresa la media (mean): "))
    params$sd <- as.numeric(readline(prompt = "Ingresa la desviación estándar (sd): "))
  } else if (dist_elegida == "binomial") {
    params$size <- as.numeric(readline(prompt = "Ingresa el número de ensayos (size): "))
    params$prob <- as.numeric(readline(prompt = "Ingresa la probabilidad de éxito (prob): "))
  } else if (dist_elegida == "poisson") {
    params$lambda <- as.numeric(readline(prompt = "Ingresa la tasa promedio (lambda): "))
  } else if (dist_elegida == "exponencial") {
    params$rate <- as.numeric(readline(prompt = "Ingresa la tasa (rate = 1/lambda): "))
  } else if (dist_elegida == "tstudent") {
    params$df <- as.numeric(readline(prompt = "Ingresa los grados de libertad (df): "))
  } else if (dist_elegida == "jicuadrado") {
    params$df <- as.numeric(readline(prompt = "Ingresa los grados de libertad (df): "))
  } else if (dist_elegida == "geometrica") {
    params$prob <- as.numeric(readline(prompt = "Ingresa la probabilidad de éxito (prob): "))
  } else if (dist_elegida == "hipergeometrica") {
    params$m <- as.numeric(readline(prompt = "Éxitos en la población (m): "))
    params$n <- as.numeric(readline(prompt = "Fracasos en la población (n): "))
    params$k <- as.numeric(readline(prompt = "Número de extracciones (k): "))
  } else if (dist_elegida == "f") {
    params$df1 <- as.numeric(readline(prompt = "Grados de libertad del numerador (df1): "))
    params$df2 <- as.numeric(readline(prompt = "Grados de libertad del denominador (df2): "))
  } else if (dist_elegida == "gamma") {
    params$shape <- as.numeric(readline(prompt = "Parámetro de forma (shape / alpha): "))
    params$rate <- as.numeric(readline(prompt = "Parámetro de tasa (rate / beta): "))
  } else if (dist_elegida == "beta") {
    params$shape1 <- as.numeric(readline(prompt = "Parámetro shape1 (alpha): "))
    params$shape2 <- as.numeric(readline(prompt = "Parámetro shape2 (beta): "))
  }
  return(params)
}

# =====================================================================
# PASO 3 Y 4: MENÚ DE OPERACIONES Y PROBABILIDADES
# =====================================================================

menu_probabilidades_avanzado <- function(dist_elegida, params, tipo_var) {
  # --- Funciones puente hacia tus fórmulas matemáticas ---
  calc_F <- function(x_val) {
    if (dist_elegida == "normal") return(mi_normal("p", x_val, params$mean, params$sd))
    if (dist_elegida == "binomial") return(mi_binomial("p", x_val, params$size, params$prob))
    if (dist_elegida == "poisson") return(mi_poisson("p", x_val, params$lambda))
    if (dist_elegida == "exponencial") return(mi_exponencial("p", x_val, params$rate))
    if (dist_elegida == "tstudent") return(mi_tstudent("p", x_val, params$df))
    if (dist_elegida == "jicuadrado") return(mi_jicuadrado("p", x_val, params$df))
    if (dist_elegida == "geometrica") return(mi_geometrico("p", x_val, params$prob))
    if (dist_elegida == "hipergeometrica") return(mi_hipergeometrico("p", x_val, params$m, params$n, params$k))
    if (dist_elegida == "f") return(mi_F("p", x_val, params$df1, params$df2))
    if (dist_elegida == "gamma") return(mi_gama("p", x_val, params$shape, params$rate))
    if (dist_elegida == "beta") return(mi_beta("p", x_val, params$shape1, params$shape2))
    return(NA)
  }
  
  calc_Q <- function(p_val) {
    if (dist_elegida == "normal") return(mi_normal("q", p_val, params$mean, params$sd))
    if (dist_elegida == "binomial") return(mi_binomial("q", p_val, params$size, params$prob))
    if (dist_elegida == "poisson") return(mi_poisson("q", p_val, params$lambda))
    if (dist_elegida == "exponencial") return(mi_exponencial("q", p_val, params$rate))
    if (dist_elegida == "tstudent") return(mi_tstudent("q", p_val, params$df))
    if (dist_elegida == "jicuadrado") return(mi_jicuadrado("q", p_val, params$df))
    if (dist_elegida == "geometrica") return(mi_geometrico("q", p_val, params$prob))
    if (dist_elegida == "hipergeometrica") return(mi_hipergeometrico("q", p_val, params$m, params$n, params$k))
    if (dist_elegida == "f") return(mi_F("q", p_val, params$df1, params$df2))
    if (dist_elegida == "gamma") return(mi_gama("q", p_val, params$shape, params$rate))
    if (dist_elegida == "beta") return(mi_beta("q", p_val, params$shape1, params$shape2))
    return(NA)
  }
  
  calc_R <- function(n_val) {
    if (dist_elegida == "normal") return(mi_normal("r", n_val, params$mean, params$sd))
    if (dist_elegida == "binomial") return(mi_binomial("r", n_val, params$size, params$prob))
    if (dist_elegida == "poisson") return(mi_poisson("r", n_val, params$lambda))
    if (dist_elegida == "exponencial") return(mi_exponencial("r", n_val, params$rate))
    if (dist_elegida == "tstudent") return(mi_tstudent("r", n_val, params$df))
    if (dist_elegida == "jicuadrado") return(mi_jicuadrado("r", n_val, params$df))
    if (dist_elegida == "geometrica") return(mi_geometrico("r", n_val, params$prob))
    if (dist_elegida == "hipergeometrica") return(mi_hipergeometrico("r", n_val, params$m, params$n, params$k))
    if (dist_elegida == "f") return(mi_F("r", n_val, params$df1, params$df2))
    if (dist_elegida == "gamma") return(mi_gama("r", n_val, params$shape, params$rate))
    if (dist_elegida == "beta") return(mi_beta("r", n_val, params$shape1, params$shape2))
    return(NA)
  }
  
  mostrar_momentos <- function() {
    cat("\n--- PROPIEDADES TEÓRICAS DEL MODELO ---\n")
    if (dist_elegida == "normal") {
      cat("Media E[X] =", params$mean, "\n")
      cat("Varianza Var(X) =", params$sd^2, "\n")
      cat("Desviación Estándar =", params$sd, "\n")
    } else if (dist_elegida == "binomial") {
      cat("Media E[X] =", params$size * params$prob, "\n")
      cat("Varianza Var(X) =", params$size * params$prob * (1 - params$prob), "\n")
      cat("Desviación Estándar =", sqrt(params$size * params$prob * (1 - params$prob)), "\n")
    } else if (dist_elegida == "poisson") {
      cat("Media E[X] =", params$lambda, "\n")
      cat("Varianza Var(X) =", params$lambda, "\n")
      cat("Desviación Estándar =", sqrt(params$lambda), "\n")
    } else if (dist_elegida == "exponencial") {
      cat("Media E[X] =", 1 / params$rate, "\n")
      cat("Varianza Var(X) =", 1 / (params$rate^2), "\n")
      cat("Desviación Estándar =", sqrt(1 / (params$rate^2)), "\n")
    } else if (dist_elegida == "tstudent") {
      cat("Media E[X] =", ifelse(params$df > 1, 0, "Indefinida (df <= 1)"), "\n")
      var_t <- ifelse(params$df > 2, params$df / (params$df - 2), ifelse(params$df > 1, Inf, "Indefinida"))
      cat("Varianza Var(X) =", var_t, "\n")
      cat("Desviación Estándar =", if(is.numeric(var_t)) sqrt(var_t) else "Indefinida", "\n")
    } else if (dist_elegida == "jicuadrado") {
      cat("Media E[X] =", params$df, "\n")
      cat("Varianza Var(X) =", 2 * params$df, "\n")
      cat("Desviación Estándar =", sqrt(2 * params$df), "\n")
    } else if (dist_elegida == "geometrica") {
      cat("Media E[X] =", (1 - params$prob) / params$prob, "\n")
      cat("Varianza Var(X) =", (1 - params$prob) / (params$prob^2), "\n")
      cat("Desviación Estándar =", sqrt((1 - params$prob) / (params$prob^2)), "\n")
    } else if (dist_elegida == "hipergeometrica") {
      N_pob <- params$m + params$n
      p_exito <- params$m / N_pob
      media_hip <- params$k * p_exito
      var_hip <- params$k * p_exito * (1 - p_exito) * ((N_pob - params$k) / (N_pob - 1))
      cat("Media E[X] =", media_hip, "\n")
      cat("Varianza Var(X) =", var_hip, "\n")
      cat("Desviación Estándar =", sqrt(var_hip), "\n")
    } else if (dist_elegida == "f") {
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
    } else if (dist_elegida == "gamma") {
      cat("Media E[X] =", params$shape / params$rate, "\n")
      cat("Varianza Var(X) =", params$shape / (params$rate^2), "\n")
      cat("Desviación Estándar =", sqrt(params$shape / (params$rate^2)), "\n")
    } else if (dist_elegida == "beta") {
      suma_ab <- params$shape1 + params$shape2
      media_beta <- params$shape1 / suma_ab
      var_beta <- (params$shape1 * params$shape2) / ((suma_ab^2) * (suma_ab + 1))
      cat("Media E[X] =", media_beta, "\n")
      cat("Varianza Var(X) =", var_beta, "\n")
      cat("Desviación Estándar =", sqrt(var_beta), "\n")
    }
  }
  
  # --- CICLO PRINCIPAL DE OPERACIONES ---
  repeat {
    cat("\n======================================================\n")
    cat(sprintf("   MENÚ DE OPERACIONES: MODELO %s\n", toupper(dist_elegida)))
    cat("======================================================\n")
    cat("1. Calcular Probabilidades (Acumulada, Supervivencia, Intervalos)\n")
    cat("2. Calcular Cuantiles (Inversa)\n")
    cat("3. Ver Valores de los Parámetros Estimados / Teóricos\n")
    cat("4. Ver Valor Esperado, Varianza y Desviación Estándar\n")
    cat("5. Generar Números Aleatorios (Simulación)\n")
    cat("0. Volver al menú anterior\n")
    
    opc_main <- readline(prompt = "Elige una opción (0-5): ")
    
    if (opc_main == "0") {
      break
    } else if (opc_main == "1") {
      # --- SUBMENÚ DE PROBABILIDADES ---
      repeat {
        cat("\n--- CÁLCULO DE PROBABILIDADES ---\n")
        cat("1. Probabilidad Acumulada: Pr(X <= a)\n")
        cat("2. Supervivencia:          Pr(X > a)\n")
        cat("3. Intervalo:              Pr(a <= X <= b)\n")
        cat("4. Valor Absoluto Menor:   Pr(|X| <= a)  [equivale a Pr(-a <= X <= a)]\n")
        cat("5. Valor Absoluto Mayor:   Pr(|X| > a)   [equivale a 1 - Pr(|X| <= a)]\n")
        cat("0. Volver al menú de operaciones\n")
        
        opc_prob <- readline(prompt = "Elige una opción (0-5): ")
        
        if (opc_prob == "0") break
        
        if (opc_prob == "1") {
          a <- as.numeric(readline(prompt = "Ingresa el valor de 'a': "))
          res <- calc_F(a)
          cat(sprintf("\nResultado -> Pr(X <= %g) = %.6f\n", a, res))
        } else if (opc_prob == "2") {
          a <- as.numeric(readline(prompt = "Ingresa el valor de 'a': "))
          res <- 1 - calc_F(a)
          cat(sprintf("\nResultado -> Pr(X > %g) = %.6f\n", a, res))
        } else if (opc_prob == "3") {
          a <- as.numeric(readline(prompt = "Ingresa el límite inferior 'a': "))
          b <- as.numeric(readline(prompt = "Ingresa el límite superior 'b': "))
          if (tipo_var == "1") res <- calc_F(b) - calc_F(a - 1)
          else res <- calc_F(b) - calc_F(a)
          cat(sprintf("\nResultado -> Pr(%g <= X <= %g) = %.6f\n", a, b, res))
        } else if (opc_prob == "4") {
          a <- as.numeric(readline(prompt = "Ingresa el valor de 'a': "))
          if (tipo_var == "1") res <- calc_F(a) - calc_F(-a - 1)
          else res <- calc_F(a) - calc_F(-a)
          cat(sprintf("\nResultado -> Pr(|X| <= %g) = %.6f\n", a, res))
        } else if (opc_prob == "5") {
          a <- as.numeric(readline(prompt = "Ingresa el valor de 'a': "))
          if (tipo_var == "1") res <- 1 - (calc_F(a) - calc_F(-a - 1))
          else res <- 1 - (calc_F(a) - calc_F(-a))
          cat(sprintf("\nResultado -> Pr(|X| > %g) = %.6f\n", a, res))
        } else {
          cat("\n❌ Opción no válida.\n")
        }
      }
    } else if (opc_main == "2") {
      p_val <- as.numeric(readline(prompt = "Ingresa la probabilidad objetivo (ej. 0.95): "))
      cat(sprintf("\nResultado -> Cuantil al %.2f%% = %g\n", p_val * 100, calc_Q(p_val)))
    } else if (opc_main == "3") {
      cat("\n--- PARÁMETROS DEL MODELO ---\n")
      for (nom in names(params)) {
        cat(sprintf("%s : %g\n", nom, params[[nom]]))
      }
    } else if (opc_main == "4") {
      mostrar_momentos()
    } else if (opc_main == "5") {
      n_val <- as.numeric(readline(prompt = "¿Cuántos números aleatorios deseas generar?: "))
      cat("\n--- NÚMEROS GENERADOS ---\n")
      print(calc_R(n_val))
    } else {
      cat("\n Opción no válida.\n")
    }
  }
}

# =====================================================================
# MOTOR COMPARATIVO VISUAL Y DE CORRELACIÓN (11 DISTRIBUCIONES)
# =====================================================================

menu_ajuste_comparativo <- function(tipo_var) {
  cat("\n--- INGRESO DE DATOS (MUESTRA) ---\n")
  cat("Puedes ingresar los números separados por comas (ej. 1.2, 3.4) \nO escribir el nombre de una variable en tu entorno (ej. vector_seca).\n")
  datos_x_str <- readline(prompt = "Datos X: ")
  
  # Lógica inteligente: Si existe una variable con ese nombre, la usa. Si no, asume que son números con comas.
  if (exists(datos_x_str)) {
    x <- get(datos_x_str)
  } else {
    x <- as.numeric(unlist(strsplit(datos_x_str, ",")))
  }
  x <- x[!is.na(x)]
  
  if (length(x) == 0) return(cat("\n❌ Error: No se detectaron datos válidos.\n"))
  
  cat("\n¿Deseas ingresar una variable Y para calcular covarianza y correlación?\n1. Sí\n2. No\n")
  if (readline(prompt = "Elige (1-2): ") == "1") {
    datos_y_str <- readline(prompt = "Datos Y: ")
    if (exists(datos_y_str)) {
      y <- get(datos_y_str)
    } else {
      y <- as.numeric(unlist(strsplit(datos_y_str, ",")))
    }
    y <- y[!is.na(y)]
    if (length(x) == length(y)) {
      cat("\n--- ESTADÍSTICA BIVARIADA ---\n")
      cat("Covarianza Cov(X,Y): ", cov(x, y), "\n")
      cat("Correlación Cor(X,Y):", cor(x, y), "\n-----------------------------\n")
    } else {
      cat("\n⚠️ Los vectores no tienen el mismo tamaño. Omitiendo correlación.\n")
    }
  }
  
  cat("\n[Ajustando modelos teóricos y generando panel en 'Plots'...]\n")
  
  m <- mean(x)
  v <- var(x)
  
  if (tipo_var == "2") { 
    # ================= CONTINUAS (7 Modelos) =================
    par(mfrow = c(2, 4)) # Cuadrícula para 7 gráficas
    x_seq <- seq(min(x), max(x), length.out = 100)
    
    # 1. Normal
    params_norm <- list(mean = m, sd = sqrt(v))
    hist(x, prob=TRUE, main="Normal", col="lightblue", xlab="X", ylab="Densidad")
    lines(x_seq, sapply(x_seq, function(val) mi_normal("d", val, params_norm$mean, params_norm$sd)), col="red", lwd=2)
    
    # 2. t de Student (Método de Momentos: Var = df / (df - 2))
    df_t <- ifelse(v > 1, max(3, round(2 * v / (v - 1))), 3) # Fallback si var <= 1
    params_t <- list(df = df_t)
    hist(x, prob=TRUE, main="t de Student", col="lightcyan", xlab="X", ylab="Densidad")
    lines(x_seq, sapply(x_seq, function(val) mi_tstudent("d", val, params_t$df)), col="red", lwd=2)
    
    # Modelos que exigen X > 0
    if (min(x) > 0) {
      # 3. Exponencial
      params_exp <- list(rate = 1 / m)
      hist(x, prob=TRUE, main="Exponencial", col="lightgreen", xlab="X", ylab="Densidad")
      lines(x_seq, sapply(x_seq, function(val) mi_exponencial("d", val, params_exp$rate)), col="red", lwd=2)
      
      # 4. Gamma
      params_gam <- list(shape = (m^2) / v, rate = m / v)
      hist(x, prob=TRUE, main="Gamma", col="lightcoral", xlab="X", ylab="Densidad")
      lines(x_seq, sapply(x_seq, function(val) mi_gama("d", val, params_gam$shape, params_gam$rate)), col="red", lwd=2)
      
      # 5. Ji-cuadrado (MoM: Media = df)
      params_chi <- list(df = max(1, round(m)))
      hist(x, prob=TRUE, main="Ji-cuadrado", col="lightgoldenrod", xlab="X", ylab="Densidad")
      lines(x_seq, sapply(x_seq, function(val) mi_jicuadrado("d", val, params_chi$df)), col="red", lwd=2)
      
      # 6. F de Snedecor (Aproximación MoM)
      df2_f <- ifelse(m > 1, max(3, round(2 * m / (m - 1))), 5)
      params_f <- list(df1 = 10, df2 = df2_f) # df1 fijado visualmente para estabilidad
      hist(x, prob=TRUE, main="F de Snedecor", col="lightpink", xlab="X", ylab="Densidad")
      lines(x_seq, sapply(x_seq, function(val) mi_F("d", val, params_f$df1, params_f$df2)), col="red", lwd=2)
      
    } else {
      cat("\n⚠️ Advertencia: Exponencial, Gamma, Ji-cuadrado y F exigen datos estrictamente positivos (X > 0). Fueron omitidas de la gráfica.\n")
    }
    
    # 7. Beta (Exige 0 < X < 1)
    if (min(x) > 0 && max(x) < 1) {
      term <- (m * (1 - m) / v) - 1
      a_beta <- max(0.1, m * term)
      b_beta <- max(0.1, (1 - m) * term)
      params_beta <- list(shape1 = a_beta, shape2 = b_beta)
      hist(x, prob=TRUE, main="Beta", col="lightyellow", xlab="X", ylab="Densidad")
      lines(x_seq, sapply(x_seq, function(val) mi_beta("d", val, params_beta$shape1, params_beta$shape2)), col="red", lwd=2)
    } else {
      plot(1, type="n", axes=FALSE, xlab="", ylab="", main="Beta (Omitida)")
      text(1, 1, "Datos fuera del\nrango (0,1)", cex=1.2, col="red")
      params_beta <- list(shape1 = 1, shape2 = 1) # Dummy params para evitar errores si el usuario la selecciona
    }
    
    par(mfrow = c(1, 1)) # Restaurar ventana
    
    # LOOP DE SELECCIÓN CONTINUA
    repeat {
      cat("\n--- SELECCIÓN DE MODELO CONTINUO ---\n")
      cat("1. Normal\n2. Exponencial\n3. t de Student\n4. Ji-cuadrado\n5. F de Snedecor\n6. Gamma\n7. Beta\n0. Salir\n")
      opc <- readline(prompt = "Elige el modelo que mejor se ajustó (0-7): ")
      
      if (opc == "1") menu_probabilidades_avanzado("normal", params_norm, tipo_var)
      else if (opc == "2") menu_probabilidades_avanzado("exponencial", params_exp, tipo_var)
      else if (opc == "3") menu_probabilidades_avanzado("tstudent", params_t, tipo_var)
      else if (opc == "4") menu_probabilidades_avanzado("jicuadrado", params_chi, tipo_var)
      else if (opc == "5") menu_probabilidades_avanzado("f", params_f, tipo_var)
      else if (opc == "6") menu_probabilidades_avanzado("gamma", params_gam, tipo_var)
      else if (opc == "7") {
        if(min(x) > 0 && max(x) < 1) menu_probabilidades_avanzado("beta", params_beta, tipo_var)
        else cat("\n❌ Los datos no están entre 0 y 1. Beta no aplicable.\n")
      }
      else if (opc == "0") break
      else cat("\n Opción no válida.\n")
    }
    
  } else { 
    # ================= DISCRETAS (4 Modelos) =================
    par(mfrow = c(2, 2)) # Cuadrícula de 2x2
    
    # Parámetros MoM / MLE
    params_poi <- list(lambda = m)
    size_bin <- max(x) 
    params_bin <- list(size = size_bin, prob = min(1, m / size_bin))
    params_geo <- list(prob = 1 / (1 + m))
    
    # Hipergeométrica (Asumimos k = max(x) y una población N razonable para visualizar)
    k_hyp <- max(x)
    N_hyp <- max(100, k_hyp * 3) # Población base
    m_hyp <- min(N_hyp, round((m * N_hyp) / k_hyp))
    params_hyp <- list(m = m_hyp, n = N_hyp - m_hyp, k = k_hyp)
    
    freq_rel <- table(x) / length(x)
    x_uni <- as.numeric(names(freq_rel))
    
    y_poi <- sapply(x_uni, function(val) mi_poisson("d", val, params_poi$lambda))
    y_bin <- sapply(x_uni, function(val) mi_binomial("d", val, params_bin$size, params_bin$prob))
    y_geo <- sapply(x_uni, function(val) mi_geometrico("d", val, params_geo$prob))
    y_hyp <- sapply(x_uni, function(val) mi_hipergeometrico("d", val, params_hyp$m, params_hyp$n, params_hyp$k))
    
    max_y <- max(c(freq_rel, y_poi, y_bin, y_geo, y_hyp), na.rm=TRUE) * 1.15 
    
    # 1. Poisson
    plot(x_uni, as.numeric(freq_rel), type="h", lwd=8, col="lightblue", ylim=c(0, max_y), main="Poisson", xlab="X", ylab="Probabilidad")
    points(x_uni, y_poi, col="red", pch=16, cex=1.5)
    # 2. Binomial
    plot(x_uni, as.numeric(freq_rel), type="h", lwd=8, col="lightgreen", ylim=c(0, max_y), main="Binomial", xlab="X", ylab="Probabilidad")
    points(x_uni, y_bin, col="red", pch=16, cex=1.5)
    # 3. Geométrica
    plot(x_uni, as.numeric(freq_rel), type="h", lwd=8, col="lightcoral", ylim=c(0, max_y), main="Geométrica", xlab="X", ylab="Probabilidad")
    points(x_uni, y_geo, col="red", pch=16, cex=1.5)
    # 4. Hipergeométrica
    plot(x_uni, as.numeric(freq_rel), type="h", lwd=8, col="lightgoldenrod", ylim=c(0, max_y), main="Hipergeométrica", xlab="X", ylab="Probabilidad")
    points(x_uni, y_hyp, col="red", pch=16, cex=1.5)
    
    par(mfrow = c(1, 1))
    
    # LOOP DE SELECCIÓN DISCRETA
    repeat {
      cat("\n--- SELECCIÓN DE MODELO DISCRETO ---\n")
      cat("1. Poisson\n2. Binomial\n3. Geométrica\n4. Hipergeométrica\n0. Salir\n")
      opc <- readline(prompt = "Elige el modelo que mejor se ajustó (0-4): ")
      
      if (opc == "1") menu_probabilidades_avanzado("poisson", params_poi, tipo_var)
      else if (opc == "2") menu_probabilidades_avanzado("binomial", params_bin, tipo_var)
      else if (opc == "3") menu_probabilidades_avanzado("geometrica", params_geo, tipo_var)
      else if (opc == "4") menu_probabilidades_avanzado("hipergeometrica", params_hyp, tipo_var)
      else if (opc == "0") break
      else cat("\n Opción no válida.\n")
    }
  }
}

# Cargar librerías necesarias
library(tidyverse)
library(jsonlite)

# 1. Definir los estados a analizar
estados_interes <- c("FL", "TX", "LA", "CO", "NV", "UT")

# 2. Función para descargar datos limitados por estado
descargar_fema_estado <- function(estado) {
  cat("Descargando datos de:", estado, "...\n")
  # Top 5,000 registros por estado para tener una muestra equilibrada y rápida
  url <- paste0("https://www.fema.gov/api/open/v2/FimaNfipClaims?$filter=state%20eq%20'", estado, "'&$top=5000")
  
  # Usamos tryCatch para evitar que el código se detenga si hay un error de conexión
  tryCatch({
    datos <- fromJSON(url)$FimaNfipClaims
    return(datos)
  }, error = function(e) {
    cat("Error al descargar el estado:", estado, "\n")
    return(NULL)
  })
}

# 3. Descargar y unir todo en una sola base de datos (map_df hace el ciclo automáticamente)
datos_brutos <- map_df(estados_interes, descargar_fema_estado)

# 4. Limpieza y preparación de la base
siniestros_limpios <- datos_brutos %>%
  select(state, dateOfLoss, amountPaidOnBuildingClaim) %>%
  drop_na() %>%
  mutate(
    fecha = as.Date(dateOfLoss),
    mes = as.numeric(format(fecha, "%m")),
    monto_pagado = as.numeric(amountPaidOnBuildingClaim),
    # Clasificación de riesgo geográfico
    perfil_riesgo = if_else(state %in% c("FL", "TX", "LA"), "Alto_Riesgo", "Bajo_Riesgo"),
    # Clasificación de riesgo temporal
    temporada = if_else(mes %in% 6:11, "Huracanes", "Seca")
  ) %>%
  filter(monto_pagado > 0) # Ignoramos reclamaciones en ceros

# 5. Pequeño resumen descriptivo para confirmar la descarga
resumen_estados <- siniestros_limpios %>%
  group_by(state, perfil_riesgo) %>%
  summarise(
    Total_Siniestros = n(),
    Monto_Promedio = mean(monto_pagado),
    Varianza = var(monto_pagado),
    .groups = "drop"
  ) %>%
  arrange(desc(Monto_Promedio))

print(resumen_estados)

# =====================================================================
# EXTRACCIÓN DE VECTORES PARA ESCENARIOS (SEVERIDAD Y FRECUENCIA)
# =====================================================================
library(dplyr)

# Asegurarnos de que los estados de alta afectación estén filtrados
datos_alta_afectacion <- siniestros_limpios %>%
  filter(state %in% c("FL", "TX", "LA"))

# ---------------------------------------------------------------------
# 1 Y 2. VECTORES CONTINUOS (SEVERIDAD: Monto pagado por siniestro)
# ---------------------------------------------------------------------

# FLORIDA
sev_fl_huracanes <- datos_alta_afectacion %>% filter(state == "FL", temporada == "Huracanes") %>% pull(monto_pagado)
sev_fl_seca      <- datos_alta_afectacion %>% filter(state == "FL", temporada == "Seca") %>% pull(monto_pagado)

# TEXAS
sev_tx_huracanes <- datos_alta_afectacion %>% filter(state == "TX", temporada == "Huracanes") %>% pull(monto_pagado)
sev_tx_seca      <- datos_alta_afectacion %>% filter(state == "TX", temporada == "Seca") %>% pull(monto_pagado)

# LOUISIANA
sev_la_huracanes <- datos_alta_afectacion %>% filter(state == "LA", temporada == "Huracanes") %>% pull(monto_pagado)
sev_la_seca      <- datos_alta_afectacion %>% filter(state == "LA", temporada == "Seca") %>% pull(monto_pagado)

# ---------------------------------------------------------------------
# 3. VECTORES DISCRETOS (FRECUENCIA: Número de siniestros diarios)
# ---------------------------------------------------------------------
# Agrupamos por fecha y estado para contar cuántos siniestros hubo cada día
frecuencia_diaria <- datos_alta_afectacion %>%
  group_by(state, fecha, temporada) %>%
  summarise(numero_siniestros = n(), .groups = "drop")

# FLORIDA (Discreto)
frec_fl_huracanes <- frecuencia_diaria %>% filter(state == "FL", temporada == "Huracanes") %>% pull(numero_siniestros)
frec_fl_seca      <- frecuencia_diaria %>% filter(state == "FL", temporada == "Seca") %>% pull(numero_siniestros)

# TEXAS (Discreto)
frec_tx_huracanes <- frecuencia_diaria %>% filter(state == "TX", temporada == "Huracanes") %>% pull(numero_siniestros)
frec_tx_seca      <- frecuencia_diaria %>% filter(state == "TX", temporada == "Seca") %>% pull(numero_siniestros)

# LOUISIANA (Discreto)
frec_la_huracanes <- frecuencia_diaria %>% filter(state == "LA", temporada == "Huracanes") %>% pull(numero_siniestros)
frec_la_seca      <- frecuencia_diaria %>% filter(state == "LA", temporada == "Seca") %>% pull(numero_siniestros)


##HACER EL DEL VALOR ESPERADO (QUIEN PIERDE MENOS O GANA MÁS)- Porcentajes (pérdida de oportunidad)



