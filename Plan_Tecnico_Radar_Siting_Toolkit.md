# Plan Técnico: Radar Siting Toolkit
### Plugin QGIS de código abierto para análisis de cobertura y emplazamiento de radares meteorológicos

**Versión:** 2.0 — reemplaza el plan técnico original
**Contexto de origen:** CREWS Caribbean 2.0 / evaluación de capacidad radar de Jamaica
**Destino:** herramienta genérica para NMHS de SIDS y estados costeros montañosos

---

## 1. Objetivo y encuadre

El toolkit responde a una pregunta operativa concreta: **dado un territorio, un conjunto de restricciones y una definición de qué área importa, ¿dónde deben ubicarse uno o más radares meteorológicos, y qué verá cada configuración?**

No es una calculadora de bloqueo. El bloqueo parcial de haz es un componente del cálculo, no el producto. El producto es evidencia geométrica reproducible que sustente una recomendación de inversión ante un gobierno o un donante.

### 1.1 Doble propósito

**Uso inmediato:** dar al experto de la misión de Jamaica la capacidad de proponer ubicaciones nuevas y defenderlas con números reproducibles.

**Uso duradero:** dejar una herramienta que otro NMHS del Caribe —Barbados, Santa Lucía, Dominica, Granada, Trinidad, Belice, Bahamas— pueda instalar desde el repositorio oficial de QGIS y aplicar a su propio territorio sin escribir código ni contratar consultores.

Estos dos propósitos entran en conflicto en el calendario. La regla que resuelve el conflicto está en la sección 9 y es innegociable.

### 1.2 Fuera de alcance

Se excluyen explícitamente, y conviene declararlo en el README para gestionar expectativas:

- Atenuación por lluvia y elección de banda de frecuencia. Es una decisión de primer orden en el trópico y **no** se resuelve con geometría.
- Propagación anómala real (ductos, gradientes medidos de refractividad). Solo se parametriza el factor de radio terrestre efectivo.
- Clutter de mar, clutter de tierra, interferencia RF, compatibilidad espectral.
- Estimación cuantitativa de precipitación, nowcasting, procesamiento de datos radar.
- Estructura de torre, radomo, carga de viento, resiliencia eléctrica.
- Costes, ciclo de vida, personal.

La herramienta cubre aproximadamente el 15% de un estudio de emplazamiento completo. El README debe decirlo así.

---

## 2. Modos de operación

### 2.1 Modo Evaluación

Entrada: capa de puntos con N sitios candidatos (o radares existentes). Salida: cobertura, bloqueo y métricas comparadas para cada uno y para la red conjunta.

Es el caso de uso cuando el experto ya tiene una lista corta.

### 2.2 Modo Descubrimiento

Entrada: máscara de terreno elegible. Salida: ráster de mérito — para cada celda candidata, un valor que indica qué tan bien vería un radar situado ahí.

Es el caso de uso cuando hay que **proponer** ubicaciones. El experto interpreta el ráster, aplica criterio propio y no cartografiable, y produce la lista corta que luego pasa por el modo Evaluación.

El modo Descubrimiento no decide. Reduce el espacio de búsqueda de un territorio entero a media docena de candidatos defendibles.

### 2.3 Modo Red

Sobre un conjunto de radares (candidatos propios más radares existentes, propios o extranjeros), calcula cobertura conjunta, redundancia y huecos residuales. Es lo que sustenta cualquier recomendación sobre número de radares.

---

## 3. Núcleo algorítmico

### 3.1 Reformulación por ángulo de horizonte

El plan original iteraba sobre ángulos de elevación discretos. Eso es innecesariamente caro y produce menos información. La reformulación:

Para cada acimut φ, se muestrea el DEM a lo largo del rayo a resolución nativa y se acumula el ángulo de horizonte máximo visto hasta cada distancia. La altura aparente del terreno respecto al radar, corregida por curvatura con radio terrestre efectivo `k·Re`:

```
Δz(s) = z(s) − h₀ − s² / (2·k·Re)
θ_terreno(s) = arctan( Δz(s) / s )
Θ(s) = np.maximum.accumulate( θ_terreno(s) )
```

donde `z(s)` es la cota del terreno a distancia de rayo `s`, y `h₀ = z(radar) + altura de torre`.

De `Θ(s)` se derivan directamente los dos productos principales:

**Elevación mínima despejada** en cada rango: `Θ(s)`, más media anchura de haz si se exige despeje total.

**Altura mínima visible** sobre el terreno local:
```
H_min(s) = h₀ + s·tan(Θ(s)) + s²/(2·k·Re) − z(s)
```

Una sola pasada `O(N)` por rayo. Sin bucle de elevaciones.

### 3.2 Bloqueo parcial de haz (lista corta)

Para los candidatos de lista corta, donde importa la fracción parcial de potencia y no solo el horizonte, se aplica la geometría de sección circular de Bech et al. (2003):

```
a(s) = s · tan(θb / 2)              # radio del haz a −3 dB
y(s) = z(s) − H_haz(s)              # terreno relativo al centro del haz

PBB = [ y·√(a²−y²) + a²·arcsin(y/a) + π·a²/2 ] / (π·a²)     para |y| ≤ a
PBB = 0                                                       para y < −a
PBB = 1                                                       para y > a

CBB = np.maximum.accumulate(PBB, axis=rango)
```

Verificación de límites: `y = 0 → PBB = 0.5`; `y = a → 1`; `y = −a → 0`. Correcto.

**El acumulado es obligatorio.** Publicar PBB por bin en lugar de CBB produce mapas erróneos en todo el sector detrás de un obstáculo. Es el error más común en implementaciones caseras.

### 3.3 Reimplementación en numpy — sin wradlib

`beam_block_frac` y `cum_beam_block_frac` son en total unas cincuenta líneas de numpy. wradlib, en cambio, arrastra scipy, matplotlib, xarray y xradar, y su propia documentación advierte que no fue diseñada como librería autocontenida, que no hay wheels para todas sus dependencias (GDAL incluido) y que la vía recomendada de instalación es conda.

Un plugin QGIS que exija conda no es instalable por un NMHS caribeño desde el gestor de plugins. La dependencia se elimina. wradlib se conserva únicamente como **referencia de validación en notebook**, fuera del plugin.

Dependencias finales: numpy y GDAL, ambas ya presentes en cualquier instalación de QGIS.

### 3.4 Proyección

Azimutal equidistante centrada en cada radar (AEQD), no UTM. A 250 km de alcance el círculo de análisis tiene 500 km de diámetro y cruza husos UTM con frecuencia, con distorsión anisótropa en los bordes. AEQD preserva distancia y acimut desde el centro exactamente, que es lo único que el cálculo necesita.

### 3.5 Muestreo del DEM

**Muestreo a resolución nativa del DEM, máximo por bin.** Con Δr = 500 m sobre un DEM de 30 m se saltan 16 celdas entre bins; una cresta afilada intermedia desaparece del análisis. Nunca interpolación bilineal en el centro del bin: la bilineal suaviza los picos, que es exactamente lo que bloquea.

**Nodata como bloqueante, o fallo ruidoso.** Si una celda sin dato se propaga como cero o como NaN silencioso, se genera visibilidad ficticia. Es el bug más probable del motor y el más difícil de detectar, porque produce resultados optimistas y visualmente plausibles. Debe haber un test explícito.

### 3.6 Rendimiento

| Operación | Volumen | Tiempo estimado |
|---|---|---|
| Un candidato, 360 acimutes, 250 km, muestreo 30 m | ~3·10⁶ muestras | 1–2 s |
| Lista corta de 6 candidatos, matriz de robustez completa | ~90 corridas | minutos |
| Barrido exhaustivo, malla de 1 km sobre Jamaica (~11.000 celdas) | ~3,3·10¹⁰ muestras | horas, ejecutable de noche |

Numba y Fortran quedan fuera. El cuello de botella es I/O de DEM y reproyección, no aritmética. La preselección por restricciones no existe por velocidad sino por realismo.

---

## 4. Figura de mérito

Esta es la decisión de diseño más consecuente del proyecto, porque es lo que un revisor va a cuestionar.

### 4.1 Métrica primaria

**Fracción del área de interés donde la altura mínima visible es ≤ H, ponderada por importancia**, evaluada en dos umbrales:

- **H = 1 km sobre el terreno** — relevante para QPE y crecidas repentinas.
- **H = 3 km sobre el terreno** — relevante para seguimiento de estructura de ciclón tropical.

Dos números por candidato, no uno. Un sitio puede ser excelente para huracanes y malo para inundaciones repentinas.

### 4.2 Penalizaciones que el agregado esconde

- **Mayor sector acimutal contiguo con CBB > 50%.** Un 20% de bloqueo repartido en pellizcos es benigno; el mismo 20% concentrado en 70° es un agujero permanente en los productos operativos.
- **Área de interés no vista por debajo de 3 km por ningún radar de la configuración.** El hueco residual que la inversión no resuelve.

### 4.3 Capa de importancia

Es parámetro de entrada obligatorio, nunca cableado. Mínimo defendible: población más cuencas de respuesta rápida. La herramienta escribe en el manifiesto qué capa se usó, porque será la primera pregunta de cualquier revisor.

**Advertencia de diseño para el usuario:** los objetivos "territorio terrestre", "aproximación marítima" y "área metropolitana" empujan hacia sitios distintos. Si no se ponderan explícitamente, la herramienta devuelve un ordenamiento que codifica una decisión de política que nadie tomó conscientemente. El README debe advertirlo.

---

## 5. Datos de elevación

### 5.1 Fuente recomendada

**Copernicus GLO-30**, 30 m, derivado de TanDEM-X (adquisiciones ~2010–2018), cobertura global. Mejor que SRTM y sin sus vacíos en terreno escarpado, que es donde más daño hacen. Descarga vía OpenTopography seleccionando un rectángulo.

### 5.2 El problema DSM/DTM y la envolvente

Lo que bloquea el haz es la superficie física real, incluidos árboles y edificios — argumento a favor del DSM. Pero GLO-30 procede de radar interferométrico, cuyo centro de fase queda *dentro* del dosel: no es DSM fiel ni DTM. Es un sesgo indeterminado que no se corrige a ciegas.

**Solución: no elegir.** Correr con GLO-30 (superficie) y con FABDEM (suelo desnudo) y reportar la envolvente. La diferencia entre ambos es la banda de incertidumbre por vegetación, cuantificada en lugar de asumida, e integrada en la matriz de robustez. FABDEM reduce el error medio absoluto en zonas forestadas de 5,15 m a 2,88 m, y en dosel denso el error mediano baja de 12,95 m a 0,45 m.

**Restricción legal:** FABDEM se distribuye bajo CC BY-NC-SA 4.0 — no comercial. Una consultoría remunerada es discutiblemente uso comercial, y la cláusula ShareAlike complica redistribuirlo con la herramienta. El dataset ofrece contacto para consultas de uso comercial. Resolver por escrito antes de usarlo en un entregable, y no incluir datos FABDEM en el repositorio bajo ninguna circunstancia.

Si no se resuelve, GLO-30 solo sigue siendo defendible; se pierde la banda de incertidumbre y se declara la limitación.

### 5.3 Descartados

ASTER GDEM (ruidoso). SRTM solo como comprobación cruzada. NASADEM y ALOS AW3D30 como terceras opiniones opcionales.

### 5.4 Techo de precisión

Ningún DEM global de 30 m resuelve el campo cercano: no ve el árbol junto a la torre ni la caseta a doscientos metros, y en los primeros kilómetros el haz va bajo y estrecho. El DEM global reduce el territorio a una lista corta; el levantamiento de horizonte en campo decide entre ellos. Ver sección 7.

### 5.5 Datum vertical

GLO-30 referencia el geoide EGM2008. Irrelevante para el bloqueo diferencial, relevante si el usuario introduce una altura de antena medida por GPS, que es elipsoidal. La herramienta registra el datum en el manifiesto y advierte si se mezclan fuentes.

---

## 6. Reproducibilidad

Para un informe que justifica inversión pública, esto vale tanto como el motor.

### 6.1 Manifiesto de ejecución

Cada corrida escribe junto a las salidas un JSON con: todos los parámetros, rutas y checksums de las capas de entrada, versión del plugin, versión de GDAL y numpy, y sello de tiempo. Media jornada de desarrollo. Cierra la objeción "¿cómo llegaste a esto?".

### 6.2 Matriz de robustez automática

El producto no es una corrida. Es un barrido:

| Dimensión | Valores |
|---|---|
| Factor de radio terrestre efectivo `k` | 1,0 / 4/3 / superrefractivo |
| Altura de torre | 10 / 15 / 20 / 25 / 30 m |
| Fuente de DEM | GLO-30 / FABDEM |

Salida: tabla de ordenamiento de candidatos por combinación. Si el orden se mantiene en las treinta combinaciones, la recomendación está blindada. Si se invierte, el equipo lo sabe antes que el revisor.

El barrido de altura de torre produce además el rendimiento marginal por metro, que es el argumento presupuestario directo ante quien firma.

### 6.3 Exportación

Tabla comparativa en CSV y figuras a resolución de publicación, escritas por el algoritmo. El análisis se va a repetir varias veces cuando cambien las suposiciones; automatizar la salida ahorra más de lo que cuesta.

---

## 7. Validación

Tres niveles, en orden de valor probatorio.

**Nivel 1 — Levantamiento de horizonte en campo.** En cada sitio de lista corta, panorámica de 360° con acimutes conocidos, o clinómetro y brújula. Comparación del perfil de horizonte observado contra el predicho por el DEM. Valida directamente la fuente de incertidumbre dominante —campo cercano, vegetación, edificaciones— y es lo único que la detecta. Horas por sitio, dentro de una misión de campo ya presupuestada.

**Nivel 2 — Contra wradlib, en notebook.** Mismo DEM, mismo sitio, comparación de CBB celda a celda. No valida la física; valida el código. Media jornada. Elimina la objeción más fácil que puede hacerse a una implementación propia.

**Nivel 3 — Sombras observadas.** Acumulación de imágenes públicas de radar durante un periodo largo: los sectores sin eco persistente revelan bloqueo real. Cualitativo (las imágenes son PNG con escala de color; los datos ya se perdieron) y no admisible como evidencia cuantitativa. Útil como figura de contraste visual si coincide.

---

## 8. Arquitectura

```
radar_siting_toolkit/
├── __init__.py
├── metadata.txt
├── plugin.py                       # registro del provider
├── LICENSE                         # GPL v2+ (obligatorio de facto en QGIS)
├── processing/
│   ├── provider.py
│   ├── alg_evaluate.py             # modo Evaluación
│   ├── alg_discover.py             # modo Descubrimiento
│   ├── alg_network.py              # modo Red
│   └── alg_robustness.py           # matriz de robustez
├── core/
│   ├── geometry.py                 # AEQD, muestreo de rayos, curvatura
│   ├── horizon.py                  # ángulo de horizonte acumulado
│   ├── blockage.py                 # Bech PBB / CBB en numpy
│   ├── metrics.py                  # figura de mérito, penalizaciones
│   ├── dem.py                      # lectura, nodata, datum, validación
│   └── manifest.py                 # trazabilidad de ejecución
├── config/
│   └── jamaica.yaml                # caso ejemplo
├── examples/
│   └── jamaica/                    # datos, config y salida reproducible
├── tests/
└── docs/                           # en inglés
```

**Sin diálogo Qt propio.** `QgsProcessingParameterFeatureSource` acepta la capa de puntos, `QgsProcessingParameterPoint` da selección por clic en el mapa, y el diálogo por lotes de Processing da la comparación multi-sitio. Una semana entera del plan original desaparece sin pérdida funcional.

**Sin descargador de DEM integrado.** Existen plugins dedicados y mantenerlo es deuda técnica. Se acepta cualquier ráster cargado y se documenta GLO-30 como fuente recomendada.

**Dirigido por configuración, no generalizado por especulación.** Un YAML por país con radares, DEM, capa de importancia y parámetros de escaneo. Jamaica se publica como caso ejemplo con datos y salida reproducible dentro del repositorio. El segundo país dirá qué hay que generalizar de verdad; anticiparlo gasta presupuesto en abstracción que nadie usa. El caso ejemplo trabajado es lo que hará que otro NMHS lo adopte.

---

## 9. Plan de trabajo

### Regla de ruta crítica

**El plugin no va en la ruta crítica del informe.** El análisis se hace en notebook, se entrega el informe, y después se convierte en plugin. Invertir el orden significa gastar en bugs de empaquetado QGIS el tiempo que se le debe a la evaluación de personal, mantenimiento y sostenibilidad, que es donde estos informes suelen fallar.

### Fases

| Fase | Contenido | Duración | Entregable |
|---|---|---|---|
| **0** | Notebook: motor de horizonte + Bech, validado contra wradlib. Análisis de Jamaica completo. | Semana 1 | Resultados para el informe. **Ruta crítica.** |
| **1** | `alg_evaluate` como QgsProcessingAlgorithm, con manifiesto. | Semana 2 | Algoritmo en la Caja de Herramientas |
| **2** | `alg_discover` con máscara de restricciones; `alg_network`. | Semana 3 | Modos de propuesta y red |
| **3** | `alg_robustness`, exportación CSV y figuras, métricas ponderadas. | Semana 4 | Producto de decisión completo |
| **4** | Empaquetado, tests, caso ejemplo Jamaica, documentación en inglés, publicación en repositorio QGIS. | Semana 5 | Plugin público |

La fase 0 alimenta el informe con independencia de si las fases 1–4 se completan a tiempo. Ese es el punto.

### Verificación de campo

El levantamiento de horizonte ocurre durante la misión, en paralelo, y realimenta la fase 0. No es una fase secuencial.

---

## 10. Riesgos

| Riesgo | Impacto | Mitigación |
|---|---|---|
| Nodata del DEM produce visibilidad ficticia | Alto — resultados erróneos plausibles | Test explícito; nodata bloqueante por defecto |
| Licencia FABDEM no resoluble | Medio — se pierde la banda de incertidumbre | GLO-30 solo, limitación declarada |
| Capa de importancia inexistente o no consensuada | Alto — invalida la figura de mérito | Definirla con MSJ antes de correr nada |
| El plugin absorbe tiempo del informe | Alto | Regla de ruta crítica, fase 0 primero |
| Objeción a implementación propia de Bech | Medio | Validación nivel 2, media jornada |
| DEM anterior a Melissa no refleja el dosel actual del occidente | Medio | Envolvente GLO-30/FABDEM; declarar en el informe |

---

## 11. Publicación

Licencia GPL v2 o posterior, obligatoria de facto para plugins QGIS. Repositorio público con el caso Jamaica reproducible. Documentación en inglés, que es el idioma operativo del Caribe anglófono. Nombre genérico — no "Jamaica" — para no señalizar que es una herramienta de un solo país.

Publicación en el repositorio oficial de plugins de QGIS al cierre de la fase 4. Un plugin que solo vive en GitHub no lo instala un NMHS.
