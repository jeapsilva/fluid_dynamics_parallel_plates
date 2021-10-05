# FORMULAÇÃO EF-LEA para HAGEN-POISSEUILLE.
# Autora: Jéssica Aparecida Silva
# Engenharia Mecânica - Universidade Federal do Rio de Janeiro

#------------------------------------------------------------------
#                                                         #
#           Bibliotecas utilizadas                        #
#                                                         #
#-------------------------------------------------------------------   

import numpy as np
import matplotlib.pyplot as plt
import meshio

#------------------------------------------------------------------
#                                                         #
#           Seção com as funções utilizadas               #
#                                                         #
#-------------------------------------------------------------------   

def get_X():
    f = open("Nodes.txt","r")
    points = []
    for i in range(0,3294):
        line = f.readline()
        line = line.strip('\n')
        line = line.split(" ")
        line = line[1:]
        for i in range(0,3):
            line[i] = float(line[i])
        points.append(line)
    f.close()
    return points

def get_IEN():
    f = open("IEN.txt","r")
    points = []
    for i in range(0,6298):
        line = f.readline()
        line = line.strip('\n')
        line = line.split(" ")
        line = line[4:]
        for i in range(0,3):
            line[i] = int(line[i]) - 1
        points.append(line)
    f.close()
    return points
  
#------------------------------------------------------------------
#                                                         #
#           Seção para leitura de malha                   #
#                                                         #
#-------------------------------------------------------------------  
    
minY = 0
maxY = 2
minX = 0
maxX = 8
    
xc = (maxX-minX)/2.0
yc = (maxY-minY)/2.0

raio = 0.50
d = 0.001

IEN = np.array(get_IEN())
points = get_X()
npoints = len(points) # numero de pontos
ne = len(IEN)   


X = np.zeros(npoints,dtype='float') 
Y = np.zeros(npoints,dtype='float') 
for i in range(0,npoints):
    X[i] = points[i][0]
    Y[i] = points[i][1]
      
# Pontos de condição de contorno
cc = np.zeros((npoints,1),dtype='float')   # pontos que possuem uma condicao de contorno de psi (Top, Bottom, Inlet)
Fc = np.zeros((npoints,1),dtype='float')   # condicoes de contorno de psi inicial (Top, Bottom, Inlet)

ccoutlet = np.zeros((npoints,1),dtype='float')   # pontos que possuem uma condicao de contorno de psi no outlet
Fcoutlet = np.zeros((npoints,1),dtype='float') 

ccw = np.zeros((npoints,1),dtype='float')   # pontos de contorno para usar em omega

cc_velo = np.zeros((npoints,1),dtype='float')   # pontos que possuem uma condicao de contorno de velocidade (Top, Bottom, Inlet)
Fu = np.zeros((npoints,1),dtype='float')        # condições de contorno de u no topo e fundo 
Fv = np.zeros((npoints,1),dtype='float')    # codições de contorno de v no topo e fundo
Fuoutlet  = np.zeros((npoints,1),dtype='float')     

# Criando as matrizes de condição de contorno de Psi e de velocidades vx e vy
for i in range(0,npoints):
                  
    # Verificando parte inferior do canal 
    if Y[i]== minY :
        cc[i] = 1.0
        Fc[i] = 0.0
        cc_velo[i] = 1.0
        Fu[i] = 0.0
        Fv[i] = 0.0
        ccw[i] = 1.0
    # Verificando parte superior do canal
    if Y[i] == maxY :
        cc[i] = 1.0
        Fc[i] = 2.0
        cc_velo[i] = 1.0
        Fu[i] = 0.0
        Fv[i] = 0.0
        ccw[i] = 1.0
    # Verificando inlet 
    if X[i] ==  minX :
        cc[i] = 1.0
        Fc[i] = Y[i]
        cc_velo[i] = 1.0
        Fu[i] = 1.0   
        Fv[i] = 0.0
        ccw[i] = 1.0
    # Verificando outlet
    if X[i] == maxX:
        ccoutlet[i] = 1.0   #cc de psi
        Fcoutlet[i] = 0.0  
        Fu[i]=0.0           #cc de velocidade
        Fv[i]=0.0
        cc_velo[i]=0.0
        ccw[i] = 1.0  #cc de omega
    # Verificando o cilindro
    if (raio -d) < ((xc-X[i])**2+ (Y[i]-yc)**2)**(0.5) < (raio + d) :
        cc[i] = 1.0
        Fc[i] = 1.0
        cc_velo[i] = 1.0
        Fu[i] = 0.0   
        Fv[i] = 0.0
#------------------------------------------------------------------
#                                                         #
#      Seção com criação das matrizes elementares         #
#                e matrizes globais                       #
#                                                         #
#------------------------------------------------------------------- 
        
# inicializando as matrizes K e M, vetor condição de contorno F 
# e os gradientes na direção X e Y

K =  np.zeros( (npoints,npoints), dtype='float')
M =  np.zeros( (npoints,npoints), dtype='float')
Gx =  np.zeros( (npoints,npoints), dtype='float')
Gy =  np.zeros( (npoints,npoints), dtype='float')

    # loop dos elementos da malha
for e in range(0,ne):
    v = IEN[e]
    det = X[v[2]]*( Y[v[0]]-Y[v[1]]) + X[v[0]]*( Y[v[1]]-Y[v[2]]) + X[v[1]]*(-Y[v[0]]+Y[v[2]])
    area = det/2.0
 
    # matriz de massa do elemento    
    m = (area/12.0) * np.array([ [2.0, 1.0, 1.0],
                                [1.0, 2.0, 1.0],
                                [1.0, 1.0, 2.0] ])
    
    b1 = Y[v[1]]-Y[v[2]]
    b2 = Y[v[2]]-Y[v[0]]
    b3 = Y[v[0]]-Y[v[1]]

    c1 = X[v[2]]-X[v[1]]
    c2 = X[v[0]]-X[v[2]]
    c3 = X[v[1]]-X[v[0]]
    
    # matriz do gradiente
    B = (1.0/(2.0*area)) * np.array([ [b1, b2, b3],
                                      [c1, c2, c3] ])
    
     #matriz do divergente
    BT = B.transpose()
    
    # matriz de rigidez do elemento
    ke = area*np.dot(BT,B)
    
    gxe = (1.0/6.0)*np.array([ [b1, b2, b3],
                               [b1, b2, b3],
                               [b1, b2, b3] ])
    gye = (1.0/6.0)*np.array([ [c1, c2, c3],
                               [c1, c2, c3],
                               [c1, c2, c3] ])
    # matrizes globais
    for i in range(0,3):
        ii = IEN[e,i]
        for j in range(0,3):
            jj = IEN[e,j]
            K[ii,jj] = K[ii,jj] + ke[i,j]
            M[ii,jj] = M[ii,jj] + m[i,j]
            Gx[ii,jj] = Gx[ii,jj] + gxe[i,j]
            Gy[ii,jj] = Gy[ii,jj] + gye[i,j]
             
# implementação do numero de iteracoes no tempo
nIter = 1000
t = 0
dt = 0.01
nu =0.001 

for j in range(0,nIter):     
    print('iteracao = ',j) # visualizar a iteração
    
    if j==0:        
        # Iniciando a velocidade        
        vx = np.zeros( (npoints,1),dtype='float')
        vy = np.zeros( (npoints,1),dtype='float')
        
        # Aplicando as condições de contorno nas velocidades
        for i in range(0,npoints):
            if cc_velo[i]==1.0:
                vx[i] = Fu[i]
                vy[i] = Fv[i]

        # MATRIZ DE DIFUSÃO ARTIFICIAL Kest
        Kest = np.zeros( (npoints,npoints), dtype='float') 
        
        for e in range(0,ne):
            vet = IEN[e]
            det = X[vet[2]]*( Y[vet[0]]-Y[vet[1]]) + X[vet[0]]*( Y[vet[1]]-Y[vet[2]]) + X[vet[1]]*(-Y[vet[0]]+Y[vet[2]])
            area = det/2.0
         
            # matrizes do elemento linear
            b1 = Y[vet[1]]-Y[vet[2]]
            b2 = Y[vet[2]]-Y[vet[0]]
            b3 = Y[vet[0]]-Y[vet[1]]
        
            c1 = X[vet[2]]-X[vet[1]]
            c2 = X[vet[0]]-X[vet[2]]
            c3 = X[vet[1]]-X[vet[0]]
            
            #   matriz de estabilização
                 #criacao da velocidade media
            v1 = IEN[e,0]
            v2 = IEN[e,1]
            v3 = IEN[e,2]
            
            vx_medio = ( vx[v1] + vx[v2] + vx[v3] )/ 3.0
            vy_medio = ( vy[v1] + vy[v2] + vy[v3] )/ 3.0
        
            kestx = ((dt/2.0)*(vx_medio/4*area))*np.array([ [vx_medio*b1*b1 + vy_medio*b1*c1 , vx_medio*b1*b2 + vy_medio*b1*c2 , vx_medio*b1*b3 + vy_medio*b1*c3],
                                                            [vx_medio*b2*b1 + vy_medio*b2*c1 , vx_medio*b2*b2 + vy_medio*b2*c2 , vx_medio*b2*b3 + vy_medio*b2*c3],
                                                            [vx_medio*b3*b1 + vy_medio*b3*c1 , vx_medio*b3*b2 + vy_medio*b3*c2 , vx_medio*b3*b3 + vy_medio*b3*c3] ])
            
            kesty = ((dt/2.0)*(vy_medio/4*area))*np.array([ [vx_medio*c1*b1 + vy_medio*c1*c1 , vx_medio*c1*b2 + vy_medio*c1*c2 , vx_medio*b1*b3 + vy_medio*c1*c3],
                                                            [vx_medio*c2*b1 + vy_medio*c2*c1 , vx_medio*c2*b2 + vy_medio*c2*c2 , vx_medio*b2*b3 + vy_medio*c2*c3],
                                                            [vx_medio*c3*b1 + vy_medio*c3*c1 , vx_medio*c3*b2 + vy_medio*c3*c2 , vx_medio*b3*b3 + vy_medio*c3*c3] ]) 
            
            #criando as matrizes GLOBAIS
            for i in range(0,3):
                ii = IEN[e,i]
                for j in range(0,3):
                    jj = IEN[e,j]
                    Kest[ii,jj] = Kest[ii,jj] + kestx[i,j] + kesty[i,j]

        # Calculando Gxvy - Gyvx
        b = np.dot(Gx,vy) - np.dot(Gy,vx)
        omega = np.linalg.solve(M,b)
        omegacc = omega.copy()   
                
        # Aplicando as condições de contorno de w 
        for i in range(0,len(cc)):
            if ccw[i] == 1.0:    
                omegacc[i] = omega[i] # so existe w na parede, no restante do dominio é nulo
                        
        # Calculo de v.\nabla\omega
        VGO = np.diagflat(vx)*Gx + np.diagflat(vy)*Gy
        
        # Iniciando solver transporte da vorticidade para omega n+1            
        LHSw = (1.0/dt)*M.copy()
        
        # Vetor do lado direito para eq. de transporte da vorticidade (omega)
        RHSw = (1.0/dt)*np.dot(M.copy(),omega) -np.dot(VGO,omega) + nu*np.dot(K.copy(),omega) #+ np.dot(Kest,omega)
        
        # imposicao das ccs para omega
        for i in range(0,len(cc)):
            if ccw[i] == 1.0:
                LHSw[i,:] = 0.0 
                LHSw[i,i] = 1.0 
                RHSw[i] = omegacc[i] 

        # SOLVER eq. transporte 
        omega = np.linalg.solve(LHSw,RHSw)
        
        # Funcao corrente K \psi = M* \omega
        RHSpsi = np.dot(M,omega)
        
        # Lado esquerdo da eq K\psi = M*\omega
        LHSpsi = K.copy()
        
        # Imposicao das ccs de \psi
        for i in range(0,len(cc)):
            if ccoutlet[i]==1.0:
                RHSpsi[i] = Fcoutlet[i]
            if cc[i] == 1.0:    
                LHSpsi[i,:] = 0.0 
                LHSpsi[i,i] = 1.0 
                RHSpsi[i] =  Fc[i] 
              
        # SOLVER da funcao corrente
        psi = np.linalg.solve(LHSpsi,RHSpsi)
        
        # encontrar as novas velocidades Mvx = Gy\psi, Mvy = -Gx\psi
        b_3 = np.dot(Gy,psi)   
        M3 = M.copy()
        vx = np.linalg.solve(M3,b_3)
        
        b_4 = np.dot(Gx,psi)
        M4 = M.copy()          
        vy = -np.linalg.solve(M4,b_4)
    
        # impor cc de velocidade nos vetores vx e vy
        for i in range(0,len(cc)):
            if cc[i] == 1.0:    
                vx[i] = Fu[i]
                vy[i] = Fv[i]
                  
        
        msh = meshio.read('novo.msh')
        xyz = msh.points
        
        IEN = msh.cells_dict['triangle']
        cells=[('triangle',IEN)]
                
        meshio.write_points_cells(
                    "./vtk/malhatres"+str(t)+".vtk",
                    xyz,
                    cells,
                    # Optionally provide extra data on points, cells, etc.
                     point_data={'u':vx,'v':vy,'psi':psi, 'omega':omega},
        #             cell_data=cell_data,
        #             field_data=field_data
                    )               
        t += 1


    if j != 0 : 
        
        # MATRIZ DE DIFUSÃO ARTIFICIAL Kest
        Kest = np.zeros( (npoints,npoints), dtype='float') 
        
        for e in range(0,ne):
            vet = IEN[e]
            det = X[vet[2]]*( Y[vet[0]]-Y[vet[1]]) + X[vet[0]]*( Y[vet[1]]-Y[vet[2]]) + X[vet[1]]*(-Y[vet[0]]+Y[vet[2]])
            area = det/2.0
         
            # matrizes do elemento linear
            b1 = Y[vet[1]]-Y[vet[2]]
            b2 = Y[vet[2]]-Y[vet[0]]
            b3 = Y[vet[0]]-Y[vet[1]]
        
            c1 = X[vet[2]]-X[vet[1]]
            c2 = X[vet[0]]-X[vet[2]]
            c3 = X[vet[1]]-X[vet[0]]
            
            #   matriz de estabilização
                 #criacao da velocidade media
            v1 = IEN[e,0]
            v2 = IEN[e,1]
            v3 = IEN[e,2]
            
            vx_medio = ( vx[v1] + vx[v2] + vx[v3] )/ 3.0
            vy_medio = ( vy[v1] + vy[v2] + vy[v3] )/ 3.0
        
            kestx = ((dt/2.0)*(vx_medio/4*area))*np.array([ [vx_medio*b1*b1 + vy_medio*b1*c1 , vx_medio*b1*b2 + vy_medio*b1*c2 , vx_medio*b1*b3 + vy_medio*b1*c3],
                                                            [vx_medio*b2*b1 + vy_medio*b2*c1 , vx_medio*b2*b2 + vy_medio*b2*c2 , vx_medio*b2*b3 + vy_medio*b2*c3],
                                                            [vx_medio*b3*b1 + vy_medio*b3*c1 , vx_medio*b3*b2 + vy_medio*b3*c2 , vx_medio*b3*b3 + vy_medio*b3*c3] ])
            
            kesty = ((dt/2.0)*(vy_medio/4*area))*np.array([ [vx_medio*c1*b1 + vy_medio*c1*c1 , vx_medio*c1*b2 + vy_medio*c1*c2 , vx_medio*b1*b3 + vy_medio*c1*c3],
                                                            [vx_medio*c2*b1 + vy_medio*c2*c1 , vx_medio*c2*b2 + vy_medio*c2*c2 , vx_medio*b2*b3 + vy_medio*c2*c3],
                                                            [vx_medio*c3*b1 + vy_medio*c3*c1 , vx_medio*c3*b2 + vy_medio*c3*c2 , vx_medio*b3*b3 + vy_medio*c3*c3] ]) 
            
            #criando as matrizes GLOBAIS
            for i in range(0,3):
                ii = IEN[e,i]
                for j in range(0,3):
                    jj = IEN[e,j]
                    Kest[ii,jj] = Kest[ii,jj] + kestx[i,j] + kesty[i,j]
        
        # Calculando Gxvy - Gyvx
        b = np.dot(Gx,vy) - np.dot(Gy,vx)
        omega = np.linalg.solve(M,b)
        omegacc = omega.copy()   
                
        #Aplicando as condições de contorno de w 
        for i in range(0,len(cc)):
            if ccw[i] == 1.0:    
                omegacc[i] = omega[i] # so existe w na parede, no restante do dominio é nulo
                        
        # Calculo de v.\nabla\omega
        VGO = np.diagflat(vx)*Gx + np.diagflat(vy)*Gy
        
        # Iniciando solver transporte da vorticidade para omega n+1            
        LHSw = (1.0/dt)*M.copy()
        
        # Vetor do lado direito para eq. de transporte da vorticidade (omega)
        RHSw = (1.0/dt)*np.dot(M.copy(),omega) -np.dot(VGO,omega) + nu*np.dot(K.copy(),omega) + np.dot(Kest,omega)
        
        # imposicao das ccs para omega
        for i in range(0,len(cc)):
            if ccw[i] == 1.0:
                LHSw[i,:] = 0.0 
                LHSw[i,i] = 1.0 
                RHSw[i] = omegacc[i] 

        # SOLVER eq. transporte 
        omega = np.linalg.solve(LHSw,RHSw)
        
        # Funcao corrente K \psi = M* \omega
        RHSpsi = np.dot(M.copy(),omega)
        
        # Lado esquerdo da eq K\psi = M*\omega
        LHSpsi = K.copy()
        
        # Imposicao das ccs de \psi
        for i in range(0,len(cc)):
            if ccoutlet[i]==1.0:
                RHSpsi[i] = Fcoutlet[i]
            if cc[i] == 1.0:    
                LHSpsi[i,:] = 0.0 
                LHSpsi[i,i] = 1.0 
                RHSpsi[i] =  Fc[i] 
              
        # SOLVER da funcao corrente
        psi = np.linalg.solve(LHSpsi,RHSpsi)
        
        # encontrar as novas velocidades Mvx = Gy\psi, Mvy = -Gx\psi
        b_3 = np.dot(Gy,psi)   
        M3 = M.copy()
        vx = np.linalg.solve(M3,b_3)
        
        b_4 = np.dot(Gx,psi)
        M4 = M.copy()            
        vy = -np.linalg.solve(M4,b_4)
    
        # impor cc de velocidade nos vetores vx e vy
        for i in range(0,len(cc)):
            if cc[i] == 1.0:    
                vx[i] = Fu[i]
                vy[i] = Fv[i]
        

        
        msh = meshio.read('novo.msh')
        xyz = msh.points
        
        IEN = msh.cells_dict['triangle']
        cells=[('triangle',IEN)]
                
        meshio.write_points_cells(
                    "./vtk/malhatres"+str(t)+".vtk",
                    xyz,
                    cells,
                    # Optionally provide extra data on points, cells, etc.
                     point_data={'u':vx,'v':vy,'psi':psi, 'omega':omega},
        #             cell_data=cell_data,
        #             field_data=field_data
                    )
                    
        t += 1