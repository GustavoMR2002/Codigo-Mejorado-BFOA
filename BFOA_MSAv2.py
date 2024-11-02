from bacteria import bacteria
from chemiotaxis import chemiotaxis
import numpy as np
import matplotlib.pyplot as plt

poblacion = []
path = "C:\\secuenciasBFOA\\multiFasta.fasta"
numeroDeBacterias = 9
numRandomBacteria = 3
iteraciones = 30
tumbo = 1                                              #numero de gaps a insertar 
nado = 3
chemio = chemiotaxis()
veryBest = bacteria(path)                #mejor bacteria   
tempBacteria = bacteria(path)            #bacteria temporal para validaciones
original = bacteria(path)                #bacteria original sin gaps
globalNFE = 0      #numero de evaluaciones de la funcion objetivo

dAttr= 0.1  # 0.1
wAttr= 0.2  # 0.2
hRep = dAttr
wRep = 10   # 10

# Para almacenar los valores de fitness en cada iteración
fitness_values = []

def clonaBest(veryBest, best):
    veryBest.matrix.seqs = np.array(best.matrix.seqs)
    veryBest.blosumScore = best.blosumScore
    veryBest.fitness = best.fitness
    veryBest.interaction = best.interaction

def validaSecuencias(path, veryBest):
    # Clona a veryBest en tempBacteria
    tempBacteria.matrix.seqs = np.array(veryBest.matrix.seqs)
    # Descartar los gaps de cada secuencia
    for i in range(len(tempBacteria.matrix.seqs)):
        tempBacteria.matrix.seqs[i] = tempBacteria.matrix.seqs[i].replace("-", "")
    
    # Valida que las secuencias originales sean iguales a las secuencias de tempBacteria
    for i in range(len(tempBacteria.matrix.seqs)):
        if tempBacteria.matrix.seqs[i] != original.matrix.seqs[i]:
            print("*****************Secuencias no coinciden********************")
            return

for i in range(numeroDeBacterias):  # Población inicial
    poblacion.append(bacteria(path))

for _ in range(iteraciones):  # Número de iteraciones
    for bacteria in poblacion:
        bacteria.tumboNado(tumbo)
        bacteria.autoEvalua()
    
    chemio.doChemioTaxis(poblacion, dAttr, wAttr, hRep, wRep)
    globalNFE += chemio.parcialNFE
    best = max(poblacion, key=lambda x: x.fitness)
    
    if (veryBest == None) or (best.fitness > veryBest.fitness):
        clonaBest(veryBest, best)
    
    # Agrega el valor de fitness de la mejor bacteria en esta iteración
    fitness_values.append(veryBest.fitness)
    
    print("interacción: ", veryBest.interaction, "fitness: ", veryBest.fitness, " NFE:", globalNFE)

    chemio.eliminarClonar(path, poblacion)
    chemio.insertRamdomBacterias(path, numRandomBacteria, poblacion)
    print("población: ", len(poblacion))

# Mostrar el genoma de la mejor bacteria al final
veryBest.showGenome()
validaSecuencias(path, veryBest)

# Graficar los resultados
plt.plot(fitness_values, label="Mejor fitness")
plt.xlabel('Iteraciones')
plt.ylabel('Fitness')
plt.title('Evolución del Fitness a lo largo de las iteraciones')
plt.legend()
plt.grid(True)
plt.show()
