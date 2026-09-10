\# Métodos de Multhopp y Vortex Lattice — Ala de envergadura finita



Implementación en \*\*MATLAB\*\* de los métodos de \*\*Multhopp\*\* y \*\*Vortex Lattice (VLM)\*\* para el análisis aerodinámico de un ala de envergadura finita.



El proyecto calcula la \*\*distribución de circulación adimensional\*\* a lo largo de la envergadura mediante dos métodos numéricos diferentes: el método de Multhopp, basado en la \*\*teoría de la línea sustentadora de Prandtl\*\*, y el método Vortex Lattice, basado en la discretización de la superficie alar mediante una malla de \*\*herraduras de torbellinos\*\*.



Los resultados obtenidos mediante ambos métodos se comparan para analizar su concordancia y estudiar la influencia de la \*\*geometría en planta\*\*, el \*\*estrechamiento\*\* y la \*\*torsión geométrica\*\* sobre la distribución de cargas aerodinámicas.



\## Metodología



El análisis se estructura en cuatro etapas:



\### 1. Definición de la geometría del ala



\- Definición de un ala simétrica de planta bitrapezoidal.

\- Envergadura total de 20 m.

\- Variación lineal de la cuerda en dos tramos.

\- Incorporación de estrechamiento mediante los parámetros `τ₁` y `τ₂`.

\- Definición de una torsión geométrica variable a lo largo de la envergadura.

\- Cálculo del ángulo de ataque efectivo de cada sección.



\### 2. Método de Multhopp



\- Formulación de la ecuación de Prandtl mediante una serie de senos.

\- Discretización de la envergadura mediante coordenadas angulares.

\- Construcción del sistema matricial de Multhopp.

\- Resolución numérica de los coeficientes de la serie.

\- Obtención de la distribución de circulación adimensional `G(y)`.

\- Análisis de la convergencia de la solución.



\### 3. Método Vortex Lattice



\- Discretización de la superficie del ala mediante una malla de paneles.

\- Representación de cada panel mediante una herradura de torbellino.

\- Colocación del segmento ligado al 25 % de la cuerda del panel.

\- Colocación del punto de control al 75 % de la cuerda.

\- Cálculo de las velocidades inducidas mediante la ley de Biot–Savart.

\- Construcción de la matriz de influencia.

\- Imposición de la condición de tangencia.

\- Resolución de las intensidades de las herraduras.

\- Obtención de la circulación total de cada sección.



\### 4. Comparación y cálculo aerodinámico



\- Comparación de las distribuciones de circulación obtenidas mediante Multhopp y VLM.

\- Análisis de las diferencias numéricas entre ambos métodos.

\- Cálculo del coeficiente de sustentación local `CL(y)`.

\- Estudio de la influencia del estrechamiento y la torsión sobre la distribución de carga.



\## Resultados



\### Distribución de circulación — Multhopp



El método de Multhopp proporciona una distribución continua de circulación mediante una serie de senos.



Debido a la torsión geométrica positiva del ala, la circulación no alcanza su máximo en el centro, sino que presenta dos máximos simétricos aproximadamente en `y = ±7 m`.



!\[Distribución de circulación mediante Multhopp](resultados/circulacion\_multhopp.png)



\### Distribución de circulación — Vortex Lattice



El método Vortex Lattice reproduce la distribución de circulación mediante una discretización bidimensional de la superficie alar.



La solución presenta igualmente dos máximos simétricos próximos a `y = ±7 m`, reproduciendo el comportamiento obtenido mediante Multhopp.



!\[Distribución de circulación mediante Vortex Lattice](resultados/circulacion\_vlm.png)



\### Comparación Multhopp — Vortex Lattice



Las distribuciones obtenidas mediante ambos métodos presentan una elevada concordancia a lo largo de toda la envergadura.



Los máximos aparecen aproximadamente en la misma posición y la diferencia entre los valores máximos de circulación obtenidos mediante ambos métodos es inferior al 3 %.



Las pequeñas diferencias observadas se deben principalmente a la naturaleza de cada método: Multhopp proporciona una solución continua mediante una serie de senos, mientras que VLM obtiene valores discretos asociados a los paneles de la malla.



!\[Comparación de la distribución de circulación](resultados/comparacion\_circulacion.png)



\### Coeficiente de sustentación local



A partir de la circulación obtenida mediante Vortex Lattice se calcula la distribución del coeficiente de sustentación local:



`CL(y) = 2·B·G(y) / c(y)`



Los valores máximos aparecen en las regiones exteriores del ala debido al efecto combinado de una circulación elevada y la reducción progresiva de la cuerda.



!\[Coeficiente de sustentación local](resultados/coeficiente\_sustentacion\_local.png)



\## Ejecución



El análisis completo se encuentra implementado en:



`src/multhopp\_vortex\_lattice.m`



El script puede ejecutarse directamente desde MATLAB y realiza de forma secuencial:



1\. Definición de los parámetros geométricos del ala.

2\. Construcción de la distribución de cuerda y torsión.

3\. Resolución del método de Multhopp.

4\. Generación de la malla Vortex Lattice.

5\. Construcción de la matriz de influencia mediante Biot–Savart.

6\. Resolución de las intensidades de las herraduras.

7\. Obtención de las distribuciones de circulación.

8\. Cálculo del coeficiente de sustentación local.

9\. Comparación de ambos métodos.

10\. Generación de las gráficas de resultados.



\## Documentación



La memoria técnica completa del proyecto, incluyendo el desarrollo teórico de ambos métodos, su implementación numérica y el análisis de resultados, está disponible en:



📄 \[Memoria técnica — Métodos de Multhopp y Vortex Lattice](docs/Memoria\_Tecnica\_Multhopp\_VLM.pdf)



\## Estructura del repositorio



```text

Multhopp-Vortex-Lattice-Method/

│

├── src/

│   └── multhopp\_vortex\_lattice.m

│

├── resultados/

│   ├── circulacion\_multhopp.png

│   ├── circulacion\_vlm.png

│   ├── comparacion\_circulacion.png

│   └── coeficiente\_sustentacion\_local.png

│

├── docs/

│   └── Memoria\_Tecnica\_Multhopp\_VLM.pdf

│

├── .gitignore

└── README.md

```



\## Tecnologías y métodos



\*\*MATLAB\*\* · \*\*Multhopp Method\*\* · \*\*Vortex Lattice Method\*\* · \*\*Lifting-Line Theory\*\* · \*\*Biot–Savart\*\* · Métodos numéricos · Aerodinámica de alas finitas



\## Autoría



Proyecto académico desarrollado en equipo en el Grado en Ingeniería Aeroespacial. Implementación y análisis numérico realizados en MATLAB.



\*\*Sandra Castresana Sánchez\*\*

