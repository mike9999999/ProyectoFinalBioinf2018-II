Primer avance del proyecto final
Propuesta del proyecto.
La propuesta para el trabajo final del curso es la de analizar datos propios obtenidos de microarreglos de expresión. Los datos son parte del trabajo de investigación del posgrado, el cual lo realiza en el Instituto nacional de medicina genómica. En concreto, el presente trabajo tendrá como objetivo analizar el efecto  agudo de la sobrenutrición en el transcriptoma del tejido adiposo de rata wistar (Rattus norvergicus). 
Materiales y métodos 
La metodología con la que se evaluará el transcriptoma consiste en la hibridación de microarreglos de expresión, en este caso se utilizó el chip Raragen 2.0 ST de Affymetrix. Este chip evalúa 28,407 transcritos de los cuales, ~22,000 son codificantes para proteínas. Cabe mencionar que en este trabajo es de interés generar un análisis de transcriptoma por contraste para determinar la expresión relativa.

Análisis bioinformático
Para realizar los diferentes análisis bioinformáticos se utilelizará el software R y paqueterías como:
Affy y Oligo : para leer los archivos tipo .CEL 
Limma: Para generar contrastes estadísticos.
Gplots: Para generar gráficos de densidad y contraste de expresión.
Ggplot2: Para generar gráficos de tipo PCA.
Ragene10sttranscriptcluster.db: Para obtener anotaciones.

 Estrategia para analizar el transcriptoma:
1.- Los datos se leerán y ordenaran por condición experimental.
2.- Se realizara un análisis de calidad y se terminará su inclusión o exclusión del estudio.
3.- Se generará una normalización de datos (raw data) para homogenizar diferencia interexperimentales.
4.- Posteriormente realizaremos los contrastes estadísticos para obtener la expresión diferencial.
5.- Una vez obtenido lo anterior, se procederá a determinar el conglomerado de transcritos diferenciales y su representación en vías metabólicas.

Cronograma
5 abril: Entrega del primer reporte que contiene la propuesta de proyecto final de la clase.
Investigación de comandos para funciones de cada paquete de R a utilizar.
26 abril: Entrega de segundo reporte que contiene el desarrollo central del proyecto con anotaciones puntuales. Realizar los análisis correspondientes a calidad y normalización así como de contrastes estadísticos.
Desarrollos de las redes de interacción metabólica.
Termino del proyecto
17 de mayo: entrega final
