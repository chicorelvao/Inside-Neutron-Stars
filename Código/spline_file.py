"""
FICHEIROS NECESSÁRIOS:
#tov.out para uma estrela
#fort.75 com a densidade de particulas

DESCRIÇÃO:

VARIÁVEIS NOS COMENTÁRIOS:
rho(b) - densidade barionica
y(i) - quantidade de particulas de cada especie

"""

from scipy import interpolate
import numpy as np 
import matplotlib.pyplot as plt


#______________________________________________________________________

#Ler o ficheiro fort.75 que contem a densidade de particulas para barions (e outros...)
# e o ficheiro tov_s1_____.out que contem a densidade barionica e o raio da estrela.
inputRho = open("/home/chico/Documents/Estagio/SplinePy/DATA/fort.75", "r") 
inputR = open("/home/chico/Documents/Estagio/SplinePy/DATA/tov_s1_m1.out", "r") 
#Transformar os valores do ficheiros em arrays
#ficheiro fort.75
fileRho = inputRho.readlines()

#ficheiro tov_s1_____.out

fileR = inputR.readlines()


#______________________________________________________________________
#selecionar a DATA útil de entrada

#ciclos para ler o documento que contem o raio
#O documento é lido por lista com linhas do tipo string.
#Depois é removido os espaço que separam as colunas e de seguida é construidas matrizes
#de numeros reais.


#file tov_s1____.out   

R = np.array([0]) #array que contem todos os raios
rhoBtov = np.array([0]) #array que contem todas as densidades barionicas do file tov_s1____.out

for line in fileR:
	lineCut = line.split(" ")

	while "" in lineCut:
		lineCut.remove("")

	R = np.vstack((R,float(lineCut[0])))
	rhoBtov = np.vstack((rhoBtov ,float(lineCut[3])))

rhoBtov = np.delete(rhoBtov, 0)
R = np.delete(R, 0)


#file fort.75
rhoBfort = np.array([0]) #array que contem todas as densidade neutroes
rhoN = np.array([0]) #array que contem todas as densidade neutroes
rhoP = np.array([0]) #array que contem todas as densidade protoes
rhoE = np.array([0])#array que contem todas as densidade eletroes
rhoM = np.array([0])#array que contem todas as densidade muoes

for line in fileRho:
	line = line.split(" ")

	while "" in line:
		line.remove("")

	rhoBfort = np.vstack((rhoBfort,float(line[0])))
	rhoN = np.vstack((rhoN,float(line[1])))
	rhoP = np.vstack((rhoP,float(line[2])))
	rhoE = np.vstack((rhoE,float(line[3])))
	rhoM = np.vstack((rhoM,float(line[4])))

rhoBfort = np.delete(rhoBfort, 0)
rhoN = np.delete(rhoN, 0)
rhoP = np.delete(rhoP, 0)
rhoE = np.delete(rhoE, 0)
rhoM = np.delete(rhoM, 0)




#______________________________________________________________________
#calculo da fração de particulas
# frac(i) = rho(i) / rho(b)

fracN = np.divide(rhoN, rhoBfort) #fração de neutroes
fracP = np.divide(rhoP, rhoBfort) #fração de protroes
fracE = np.divide(rhoE, rhoBfort) #fração de eletroes
fracM = np.divide(rhoM, rhoBfort) #fração de muoes



#______________________________________________________________________



realX = R #definição do eixo do x que queremos: Raio
resolX = 0.01 #resolução do x
midY = rhoBtov # definição do valor que faz a ligação entre o raio e a fração de particulas: a densidade barionica do tov_s1____.out
midX = rhoBfort # o mesmo que a variavel midY. Desta vez, é a densB de fort.75. A relação de r e frac, define se entre as variaveis rhoBtov e rhoBfort

#spline entre o raio e rho(b)tov, return dos coeficientes nos respectivos intervalos a aplicar
coef_Intp = interpolate.splrep(realX, midY, s=0)

#construção de um eixo do x relativamente contínuo. O ultimo parametro define a resolução

contX = np.arange(R[0], R[-1], resolX)
#interpolação para todos os valores do eixo do x referido e do mid Y. Ou seja, construção do eixo do x para a proxima interpolação
# x = R / y = rhoBtov = middleX
middleX = interpolate.splev(contX, coef_Intp, der=0)


#selecionador da fraçao atual da liprint(realY)sta fracs
saveY = [] #lista que guarda dos realY
fracs = [fracN, fracP, fracE, fracM] 
for currFrac in fracs:
	
	realY = currFrac#definição do eixo do y que queremos: fração de particulas


	#spline entre o rhoBfort e frac, return dos coeficientes nos respectivos intervalos a aplicar
	coef_Sec_Intp= interpolate.splrep(midX, realY, s=0) #certo

	#interpolação para todos os valores do eixo de um novo eixo do x. Ou seja, construção de um novo eixo do y
	contY= interpolate.splev(middleX, coef_Sec_Intp, der=0)
	saveY.append(contY)
	



outputFile = open("/home/chico/Documents/Estagio/SplinePy/out.txt", "w")

for step in range(0,saveY[0].size):
	outputFile.write(str(saveY[0][step]) + "  " + str(saveY[1][step]) + "  "+ str(saveY[2][step]) + "  " + str(saveY[3][step]) + "\n")

outputFile.close()

#construção do gráfico

xLine = contX

#valores finais
#x = r / y = y(i)
plt.plot(xLine, saveY[0][::-1], label='Neutroes')
plt.plot(xLine, saveY[1][::-1], label='Protroes')
plt.plot(xLine, saveY[2][::-1], label='eletroes')
plt.plot(xLine, saveY[3][::-1], label='muoes')
legend = plt.legend(loc='best', shadow=False, fontsize='x-large')

plt.show()

#-------------------muito importante---------------
#fechar o file
inputRho.close()
