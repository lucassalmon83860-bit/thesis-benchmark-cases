# -*- coding: utf-8 -*-
"""
Auteur: Lucas Salmon
Date: 20/03/24    

Objectif : ce script insère de manière automatisée des zones cohésives d'épaiseur nulles dans un maillage Abaqus (CPE3T ou CPS3T).

Pour cela les étapes suivantes sont réalisées :
- Lire le jdd Abaqus (fichier .inp) afin de récupérer la table des noeuds et la table de connexion. 
- Dans la zone (rectangle) définie par l'utilisateur, les noeuds sont dédoublés de telle sorte que chaque élément soit définit par 3 
noeuds uniques (aucun noeud n'est utilisé dans 2 éléments). Une fois les noeuds dédoublés, la table de noeud et la table de connexion
sont reconstruites.
- Pour construire les éléments cohésifs, on isole les 4 noeuds de chaque surface de contact élément triangulare/élément triangulaire. 
Ensuite une manipulation expliquée dans le doc verification_sens_horaire_pour_element_CZM.pdf est réalisée afin de définir chaque 
élément cohésif en suivant la bonne rotation (sinon l'aire devient négative lors du chargement en traction) puis on complète la table 
de connexion avec les éléments cohésifs créés.
- Si l'utilisateur souhaite réaliser un calcul thermo-nécanique, le script créé les surfaces et intéraction de surfaces associées à chaque 
élément cohésif.
- Réécriture du JDD avec tous l'éléments nécessaires associés à l'insertion des éléments cohésifs.
"""

import numpy as np
import matplotlib.pyplot as plt
import random as rdm
import math

print("-----------------------------------------")
print("                  Début")
print("-----------------------------------------\n")

############################################### Input du script ###############################################
file_name = "plaque_entaille.inp"
output_file_name ="plaque_entaillee_CZM.inp"

# On définit un domain restrein dans lequel on souhaite insérer les zones cohésives 
"""-> attention si on ne souhaite pas intégrer des éléments cohésifs sur le domaine entier il y a un risque que 
les noeuds proches des bords du sous-domaine soient dédoublés pour rien"""
borne_x = [-101,101]
borne_y = [-101,101]

# On renseignes les données nécessaires pour la construction du matériau associé aux zones cohésives 
Depvar = 14  # nombre de variables d'états (utilisées pour stocker différents paramètres d'intéret des zones cohésives et les visualiser en post-traitement) 
liste_parametres = [438394064., 438394064., 184586974., 0., 1.e-8, 100, 0.35, 200., 45000., 2., 0.1]
gc_weibull = 'no' # 'no'  :  création d'un métériau unique pour les éléments cohésifs
                   # 'yes' :  création d'un matériau pour chaque élément cohésif avec une distribution de weibull sur gc
cv = 0.02  # coefficient de variation de la distribution de Weibull (utilisé si gc_weibull = 'yes')

# On précise si on fait un calcul mécanique ou thermo-mécanique
""" -> si on fait un calcul mécanique on n'a pas besoin de définir la gestion de l'échange thermique à travers 
chaque zone cohésive sinon on doit créer les différentes surfaces d'échange ainsi que l'interaction entre ces 
surfaces d'échange"""
simulation = 'meca' # 'meca'        : calcul mécanique
                    # 'thermo-meca' : calcul thermo-mécanique


######################################### Lecture du fichier d'entrée ##########################################
table = []
table_noeuds_pastille = [] # Création d'une table vide
table_connexion_pastille = [] # Création d'une table vide
ind = 0

with open(file_name,"r",encoding='utf-8') as obj_fichier :
    while 1 :
        ligne = obj_fichier.readline()
        if not(ligne) :
            break
        else:
            ligne1=ligne.replace('\n','')
            ligne2=ligne1.replace(' ','')
            elements = ligne2.split(',')

            # On recherche les mots clefs d'intéret
            if elements[0] == '*Node' :
              ind = 1
            elif elements[0] == '*Element' :
              ind = 2
            elif elements[0] == '*Nset' :
              ind = 0
            elif elements[0] == '*EndPart' :
              ind = 0

            # On récupère les données de la ligne si elle nous intéresse (noeuds ou éléments)
            if ind == 0 or elements[0] == '*Node' or elements[0] == '*Element' :
              table.append(ligne)
            elif ind == 1 and elements[0]!= '*Node' :
              table_vide = []
              for i in elements :
                table_vide.append(float(i))
              table_vide[0] = int(table_vide[0])
              table_noeuds_pastille.append(table_vide)
            elif ind == 2 and elements[0]!= '*Element' :
              table_vide = []
              for i in elements :
                table_vide.append(int(i))
              table_connexion_pastille.append(table_vide)

obj_fichier.close()
print(" - Lecture du fichier -> ok")

########################### Dédoublement des noeuds et reconstruction des tables ##########################
new_table_noeuds_pastille = table_noeuds_pastille
new_table_connexion_pastille = table_connexion_pastille
indice_new_noeud = len(table_noeuds_pastille)+1
table_used_nodes = []
table_duplicated_nodes = []
for i in range (len(table_connexion_pastille)):
  for j in range (1,len(table_connexion_pastille[0])):
    node = table_connexion_pastille[i][j]
    x_node = table_noeuds_pastille[node-1][1]
    y_node = table_noeuds_pastille[node-1][2]
    
    #Restriction de la zone sur laquelle on dédouble les noeuds pour insérer les zones cohésives
    if x_node>borne_x[0] and x_node<borne_x[1] and y_node>borne_y[0] and y_node<borne_y[1] :
        # on regarde si le noeud est déjà apparu dans un élément précédent
        if node in table_used_nodes :
          # on recherche l'indice du noeud dans la table 'table_noeuds_pastille'
          for k in range(len(table_noeuds_pastille)):
            if node == table_noeuds_pastille[k][0]:
              node_index = k
              break
    
          # on duplique le noeud
          new_table_noeuds_pastille.append([indice_new_noeud,table_noeuds_pastille[node_index][1],table_noeuds_pastille[node_index][2]])
          new_table_connexion_pastille[i][j] = indice_new_noeud
          indice_new_noeud += 1
    
          # on associe le nouveau noeud créé à l'ancien dans une nouvelle table
          for k in range(len(table_duplicated_nodes)):
            if node == table_duplicated_nodes[k][0]:
              node_index = k
    
              break
          table_duplicated_nodes[node_index].append(indice_new_noeud-1)
        else :
          table_used_nodes.append(node)
          table_duplicated_nodes.append([node])
        table_used_nodes.append(table_connexion_pastille[i][j])

nb_elem = len(new_table_connexion_pastille)
print(" - Dédoublement des noeuds -> ok")

########################### Création des éléments cohésifs et complétion de la table de connexion ##########################
table_couples = []
table_couples_elements = []
for i in table_duplicated_nodes :
  for j in i :
    for k in i[i.index(j)+1:]:
      # On cherche les éléments dans lesquelles ont retrouve les noeuds dédoublés
      stop = 0
      for m in new_table_connexion_pastille :
        if j in m[1:len(m)] :
          elem_j = m[0]
          stop+= 1
        elif k in m[1:len(m)] :
          elem_k = m[0]
          stop+=1
        if stop ==2 :
          table_couples.append([j,k])
          table_couples_elements.append([elem_j,elem_k])
          break
print(" - Création des tables pour construire les éléments coésifs -> ok")

# On associe chaque couple entre eux
table_couple_element_new = []
table_surface_interaction = []
increment = len(new_table_connexion_pastille)+1
compteur = 1
homothetie = 0.05
for i in table_couples_elements[:len(table_couples_elements)-2]: 
  for j in table_couples_elements[compteur:] :
    if j == i :
      ind1 = table_couples_elements.index(i)
      ind2 = ind1 + table_couples_elements[table_couples_elements.index(i)+1:].index(j) + 1
      
      # on verifie le sens de rotation de l'élément (anti-horaire necessaire)
      if len(new_table_connexion_pastille[0]) == 4 :
          elem1_node1_index = new_table_connexion_pastille[i[0]-1][1]
          elem1_node2_index = new_table_connexion_pastille[i[0]-1][2]
          elem1_node3_index = new_table_connexion_pastille[i[0]-1][3]
          x_elem1_center = new_table_noeuds_pastille[elem1_node1_index-1][1]+(2/3)*(((new_table_noeuds_pastille[elem1_node2_index-1][1]+new_table_noeuds_pastille[elem1_node3_index-1][1])/2)-new_table_noeuds_pastille[elem1_node1_index-1][1])
          y_elem1_center = new_table_noeuds_pastille[elem1_node1_index-1][2]+(2/3)*(((new_table_noeuds_pastille[elem1_node2_index-1][2]+new_table_noeuds_pastille[elem1_node3_index-1][2])/2)-new_table_noeuds_pastille[elem1_node1_index-1][2])
          
          elem2_node1_index = new_table_connexion_pastille[i[1]-1][1]
          elem2_node2_index = new_table_connexion_pastille[i[1]-1][2]
          elem2_node3_index = new_table_connexion_pastille[i[1]-1][3]
          x_elem2_center = new_table_noeuds_pastille[elem2_node1_index-1][1]+(2/3)*(((new_table_noeuds_pastille[elem2_node2_index-1][1]+new_table_noeuds_pastille[elem2_node3_index-1][1])/2)-new_table_noeuds_pastille[elem2_node1_index-1][1])
          y_elem2_center = new_table_noeuds_pastille[elem2_node1_index-1][2]+(2/3)*(((new_table_noeuds_pastille[elem2_node2_index-1][2]+new_table_noeuds_pastille[elem2_node3_index-1][2])/2)-new_table_noeuds_pastille[elem2_node1_index-1][2])
          
          CZM_elem1_node1 = [new_table_noeuds_pastille[table_couples[ind2][0]-1][1]+(homothetie-1)*(new_table_noeuds_pastille[table_couples[ind2][0]-1][1]-x_elem1_center),new_table_noeuds_pastille[table_couples[ind2][0]-1][2]+(homothetie-1)*(new_table_noeuds_pastille[table_couples[ind2][0]-1][2]-y_elem1_center)]
          CZM_elem1_node2 = [new_table_noeuds_pastille[table_couples[ind1][0]-1][1]+(homothetie-1)*(new_table_noeuds_pastille[table_couples[ind1][0]-1][1]-x_elem1_center),new_table_noeuds_pastille[table_couples[ind1][0]-1][2]+(homothetie-1)*(new_table_noeuds_pastille[table_couples[ind1][0]-1][2]-y_elem1_center)]
          CZM_elem2_node1 = [new_table_noeuds_pastille[table_couples[ind1][1]-1][1]+(homothetie-1)*(new_table_noeuds_pastille[table_couples[ind1][1]-1][1]-x_elem2_center),new_table_noeuds_pastille[table_couples[ind1][1]-1][2]+(homothetie-1)*(new_table_noeuds_pastille[table_couples[ind1][1]-1][2]-y_elem2_center)]
          CZM_elem2_node2 = [new_table_noeuds_pastille[table_couples[ind2][1]-1][1]+(homothetie-1)*(new_table_noeuds_pastille[table_couples[ind2][1]-1][1]-x_elem2_center),new_table_noeuds_pastille[table_couples[ind2][1]-1][2]+(homothetie-1)*(new_table_noeuds_pastille[table_couples[ind2][1]-1][2]-y_elem2_center)]
          
          x_CZM_center = (CZM_elem1_node1[0]+CZM_elem2_node1[0])/2
          y_CZM_center = (CZM_elem1_node1[1]+CZM_elem2_node1[1])/2
          
      elif len(new_table_connexion_pastille[0]) == 5 :
          elem1_node1_index = new_table_connexion_pastille[i[0]-1][1]
          elem1_node3_index = new_table_connexion_pastille[i[0]-1][3]
          x_elem1_center = (new_table_noeuds_pastille[elem1_node1_index-1][1]+new_table_noeuds_pastille[elem1_node3_index-1][1])/2
          y_elem1_center = (new_table_noeuds_pastille[elem1_node1_index-1][2]+new_table_noeuds_pastille[elem1_node3_index-1][2])/2
          
          elem2_node1_index = new_table_connexion_pastille[i[1]-1][1]
          elem2_node3_index = new_table_connexion_pastille[i[1]-1][3]
          x_elem2_center = (new_table_noeuds_pastille[elem2_node1_index-1][1]+new_table_noeuds_pastille[elem2_node3_index-1][1])/2
          y_elem2_center = (new_table_noeuds_pastille[elem2_node1_index-1][2]+new_table_noeuds_pastille[elem2_node3_index-1][2])/2
          
      CZM_elem1_node1 = [new_table_noeuds_pastille[table_couples[ind2][0]-1][1]+(homothetie-1)*(new_table_noeuds_pastille[table_couples[ind2][0]-1][1]-x_elem1_center),new_table_noeuds_pastille[table_couples[ind2][0]-1][2]+(homothetie-1)*(new_table_noeuds_pastille[table_couples[ind2][0]-1][2]-y_elem1_center)]
      CZM_elem1_node2 = [new_table_noeuds_pastille[table_couples[ind1][0]-1][1]+(homothetie-1)*(new_table_noeuds_pastille[table_couples[ind1][0]-1][1]-x_elem1_center),new_table_noeuds_pastille[table_couples[ind1][0]-1][2]+(homothetie-1)*(new_table_noeuds_pastille[table_couples[ind1][0]-1][2]-y_elem1_center)]
      CZM_elem2_node1 = [new_table_noeuds_pastille[table_couples[ind1][1]-1][1]+(homothetie-1)*(new_table_noeuds_pastille[table_couples[ind1][1]-1][1]-x_elem2_center),new_table_noeuds_pastille[table_couples[ind1][1]-1][2]+(homothetie-1)*(new_table_noeuds_pastille[table_couples[ind1][1]-1][2]-y_elem2_center)]
      CZM_elem2_node2 = [new_table_noeuds_pastille[table_couples[ind2][1]-1][1]+(homothetie-1)*(new_table_noeuds_pastille[table_couples[ind2][1]-1][1]-x_elem2_center),new_table_noeuds_pastille[table_couples[ind2][1]-1][2]+(homothetie-1)*(new_table_noeuds_pastille[table_couples[ind2][1]-1][2]-y_elem2_center)]
          
      x_CZM_center = (CZM_elem1_node1[0]+CZM_elem2_node1[0])/2
      y_CZM_center = (CZM_elem1_node1[1]+CZM_elem2_node1[1])/2
          
      CZM_elem1_node1_new = [CZM_elem1_node1[0]-x_CZM_center, CZM_elem1_node1[1]-y_CZM_center]
      CZM_elem1_node2_new = [CZM_elem1_node2[0]-x_CZM_center, CZM_elem1_node2[1]-y_CZM_center]
          
      if CZM_elem1_node1_new[1]*CZM_elem1_node2_new[1]>0 :
          point_fictif = [-10.,0.]
      else :
          point_fictif = [0.,-10.]
          if CZM_elem1_node1_new[0]*CZM_elem1_node2_new[0]<0 :
             print('Problème avec le choix du point fictif ! ')
             print(CZM_elem1_node1_new[0],CZM_elem1_node2_new[0])
          
      distance_fictive_carre = point_fictif[0]**2+point_fictif[1]**2
      distance1_carre = CZM_elem1_node1_new[0]**2 + CZM_elem1_node1_new[1]**2
      distance1bis_carre = (CZM_elem1_node1_new[0]-point_fictif[0])**2 + (CZM_elem1_node1_new[1]-point_fictif[1])**2
      distance2_carre = CZM_elem1_node2_new[0]**2 + CZM_elem1_node2_new[1]**2
      distance2bis_carre = (CZM_elem1_node2_new[0]-point_fictif[0])**2 + (CZM_elem1_node2_new[1]-point_fictif[1])**2

      angle1 = np.arccos((distance_fictive_carre + distance1_carre - distance1bis_carre)/(2.*np.sqrt(distance_fictive_carre)*np.sqrt(distance1_carre)))
      angle2 = np.arccos((distance_fictive_carre + distance2_carre - distance2bis_carre)/(2.*np.sqrt(distance_fictive_carre)*np.sqrt(distance2_carre)))
      if CZM_elem1_node1_new[1]*CZM_elem1_node2_new[1]>0 :
          if CZM_elem1_node1_new[1] > 0 :
              angle = angle1 - angle2
          else :
              angle = angle2 - angle1
      else :
          if CZM_elem1_node1_new[0] < 0 :
              angle = angle1 - angle2
          else :
              angle = angle2 - angle1
  

      if angle > 0 :
          sens = 1
      else :
          sens = -1

      if sens == 1 :
          element = [increment,table_couples[ind2][0],table_couples[ind1][0],table_couples[ind1][1],table_couples[ind2][1]]
      else :
          element = [increment,table_couples[ind1][0],table_couples[ind2][0],table_couples[ind2][1],table_couples[ind1][1]]
      new_table_connexion_pastille.append(element)
      table_couple_element_new.append(j)
      table_surface_interaction.append([[table_couples[ind1][0],table_couples[ind2][0]],[table_couples[ind1][1],table_couples[ind2][1]]])
      increment+=1
  compteur+=1
nb_elem_czm = len(new_table_connexion_pastille)- nb_elem
print(" - Table de connexion des éléments cohésifs -> ok")

########################### On créé les Elset et Surfaces pour interaction ##########################
if simulation == 'meca' : 
  print(" - Création des Elset et Surfaces -> pas besoin pour une simulation mécanique")
elif simulation =='thermo-meca':
  ind = 1
  ind_elset = 1
  ind_surface = 1
  table_Elset_interaction = []
  table_Surfaces_interaction = []
  liste_Elsetij_m = []
  for i in table_surface_interaction :
    i_index = table_surface_interaction.index(i)
    
    for j in i :
      j_index = i.index(j)
      # matrice
      Elsetij_m = table_couple_element_new[i_index][j_index]
      position_node1 = new_table_connexion_pastille[Elsetij_m-1][1:].index(j[0])
      position_node2 = new_table_connexion_pastille[Elsetij_m-1][1:].index(j[1])
      order_nodes = [position_node1,position_node2]
      if len(new_table_connexion_pastille[0]) == 4 :
          if order_nodes == [0,1] or order_nodes == [1,0]:
            surface = 'S1'
          elif order_nodes == [1,2] or order_nodes == [2,1]:
            surface = 'S2'
          elif order_nodes == [2,0] or order_nodes == [0,2]:
            surface = 'S3'
      elif len(new_table_connexion_pastille[0]) == 5 :
          if order_nodes == [0,1] or order_nodes == [1,0]:
            surface = 'S1'
          elif order_nodes == [1,2] or order_nodes == [2,1]:
            surface = 'S2'
          elif order_nodes == [2,3] or order_nodes == [3,2]:
            surface = 'S3'
          elif order_nodes == [3,0] or order_nodes == [0,3]:
            surface = 'S4'
      
      new_elset = ['m_set-'+str(ind_elset),Elsetij_m]
      if Elsetij_m in liste_Elsetij_m :
          nb_elset = liste_Elsetij_m.index(Elsetij_m)
          table_Surfaces_interaction.append(['m_surf-'+str(ind_surface)+'-'+str(surface),table_Elset_interaction[nb_elset][0],surface])
      else :
          table_Elset_interaction.append(new_elset)
          liste_Elsetij_m.append(Elsetij_m)
          table_Surfaces_interaction.append(['m_surf-'+str(ind_surface)+'-'+str(surface),new_elset[0],surface])
          ind_elset+=1
      
      ind_surface+=1
  print(" - Création des Elset et Surfaces -> ok")

########################### On réécrit le fichier d'entrée Abaqus ##########################
ind = 'null'
module = 'null'
new2_table_noeuds_pastille =new_table_noeuds_pastille
with open(output_file_name,"w",encoding='utf-8') as obj_fichier :
  for i in table:

    if ind == 'Node':
      for j in new2_table_noeuds_pastille :
        line = ' '+str(j[0])+', '+str(j[1])+', '+str(j[2])+'\n'
        obj_fichier.write(line)
        ind = 'null'

    elif ind == 'Element' :
      for j in new_table_connexion_pastille[:nb_elem] :
        line = ' '
        for k in j[:len(j)-1] :
            line += str(k)+', '
        line += str(j[len(j)-1])+'\n'
        obj_fichier.write(line)
      # ajout des elements czm
      obj_fichier.write('*Element, type=COH2D4T\n')
      for j in new_table_connexion_pastille[nb_elem:] :
        line = ' '+str(j[0])+', '+str(j[1])+', '+str(j[2])+', '+str(j[3])+', '+str(j[4])+'\n'
        obj_fichier.write(line)
      ind = 'null'

    elif ind == 'Nset' :
      if module == 'PARTS' :
          line = ' 1, '+str(len(new_table_noeuds_pastille))+', 1\n'
          obj_fichier.write(line)
      if module == 'ASSEMBLY' :
          # on rajoute les noeuds dédoublés
          ligne=i.replace('\n','')
          ligne1=ligne.replace(' ','')
          elements = ligne1.split(',')
          new_ligne = ligne
          if len(elements) > 2 :
              for j in elements :
                  for k in table_duplicated_nodes :
                      if int(j) == k[0] :
                          if len(k) > 1 :
                              for m in k[1:] :
                                  new_ligne+= ', '+str(m)
          new_ligne+='\n'
          obj_fichier.write(new_ligne)
          obj_fichier.write('**\n') 
    
    elif ind == 'Elset' :
      if module == 'PARTS' :
          line = ' 1, '+str(nb_elem)+', 1\n'
          obj_fichier.write(line)

    elif ind == 'Solid Section' :
      obj_fichier.write(',\n')
      if gc_weibull == 'no' :
        obj_fichier.write('*Nset, nset=Set-2, generate\n')
        line = ' 1, '+str(len(new_table_noeuds_pastille))+', 1\n'
        obj_fichier.write(line)

        obj_fichier.write('*Elset, elset=Set-2, generate\n')
        line = ' '+str(nb_elem+1)+', '+str(len(table_connexion_pastille))+', 1\n'
        obj_fichier.write(line)

        obj_fichier.write('** Section: Section-2\n')
        obj_fichier.write('*Cohesive Section, elset=Set-2, material=Material-2, response=TRACTION SEPARATION\n')

      elif gc_weibull == 'yes' :  
        for j in new_table_connexion_pastille[nb_elem:]:
            line = ' '+str(j[1])+', '+str(j[2])+', '+str(j[3])+', '+str(j[4])+'\n'
            obj_fichier.write('*Nset, nset=Set-CZM'+str(j[0])+'\n')
            obj_fichier.write(line)
      
            obj_fichier.write('*Elset, elset=Set-CZM'+str(j[0])+'\n')
            line = ' '+str(j[0])+'\n'
            obj_fichier.write(line)
      
            obj_fichier.write('** Section: Section-CZM'+str(j[0])+'\n')
            obj_fichier.write('*Cohesive Section, elset=Set-CZM'+str(j[0])+', material=Material-CZM'+str(j[0])+', response=TRACTION SEPARATION\n')
            obj_fichier.write(',\n')
      ind = 'null'
    
    elif ind == 'Set' :
      obj_fichier.write(i)
      
      # initialisation du Set-3 qui regroupe tous le domaine
      obj_fichier.write('*Nset, nset=Set-3, instance=Part-1-1, generate\n')
      line = ' 1, '+str(len(new_table_noeuds_pastille))+', 1\n'
      obj_fichier.write(line)
      obj_fichier.write('*Elset, elset=Set-3, instance=Part-1-1, generate\n')
      line = ' 1, '+str(len(table_connexion_pastille))+', 1\n'
      obj_fichier.write(line)
      obj_fichier.write('**\n') 
      
      if simulation == 'thermo-meca' :
        # on construit tous les Elset ( théoriquement un par élémént TRI3 mais attention doublons pour faciliter la construction du fichier!)
        for j in table_Elset_interaction : 
          obj_fichier.write('*Elset, elset='+j[0]+', internal, instance=Part-1-1, generate\n')
          obj_fichier.write(' '+str(j[1])+', '+str(j[1])+', 1\n')
        obj_fichier.write('**\n')
        
        # on construit toutes les surfaces pour les intéractions
        for j in table_Surfaces_interaction :
            obj_fichier.write('*Surface, type=ELEMENT, name='+j[0]+'\n')
            index_j= table_Surfaces_interaction.index(j)
            obj_fichier.write(j[1]+', '+j[2]+'\n')
        obj_fichier.write('**\n')
      
    elif ind == 'Material2' :
      obj_fichier.write(i)
      # définition du materiau CZM si on ne souhaite pas de distribution de Weibull sur gc
      if gc_weibull == 'no' : 
        obj_fichier.write('*Material, name=Material-2\n')
        obj_fichier.write('*Depvar\n')
        obj_fichier.write(' '+str(Depvar)+',\n')
        obj_fichier.write('*User Material, constants='+str(len(liste_parametres))+'\n')
        obj_fichier.write(str(liste_parametres[0])+', '+str(liste_parametres[1])+', '+str(liste_parametres[2])+', '+str(liste_parametres[3])+', '+str(liste_parametres[4])+', '+str(liste_parametres[5])+', '+str(liste_parametres[6])+', '+str(liste_parametres[7])+'\n')
        obj_fichier.write(str(liste_parametres[8])+', '+str(liste_parametres[9])+', '+str(liste_parametres[10])+'\n')
        obj_fichier.write('**\n') 
      # définition du materiau CZM si on souhaite une distribution de Weibull sur gc
      elif gc_weibull == 'yes' :
        m = 2
        Gc_moy = liste_parametres[6]
        Gc_0 = 2.1586552*Gc_moy*cv
        Gc_min = Gc_moy*(1-(1.9130584*cv))
        for j in new_table_connexion_pastille[nb_elem:]:
            rand = rdm.random()
            Gc = round(((-math.log(rand))**(1/m))*Gc_0+Gc_min,4)
            obj_fichier.write('*Material, name=Material-CZM'+str(j[0])+'\n')
            obj_fichier.write('*Depvar\n')
            obj_fichier.write(' '+str(Depvar)+',\n')
            obj_fichier.write('*User Material, constants='+str(len(liste_parametres))+'\n')
            obj_fichier.write(str(liste_parametres[0])+', '+str(liste_parametres[1])+', '+str(liste_parametres[2])+', '+str(liste_parametres[3])+', '+str(liste_parametres[4])+', '+str(liste_parametres[5])+', '+str(Gc)+', '+str(liste_parametres[7])+'\n')
            obj_fichier.write(str(liste_parametres[8])+', '+str(liste_parametres[9])+', '+str(liste_parametres[10])+'\n')
            obj_fichier.write('**\n') 
      
      # définition des propriétés des intéractions
      obj_fichier.write('** INTERACTION PROPERTIES\n') 
      obj_fichier.write('**\n') 
      obj_fichier.write('*Surface Interaction, name=IntProp-1\n') 
      obj_fichier.write('1.,\n') 
      obj_fichier.write('*Gap Conductance\n') 
      obj_fichier.write(' 3e6,  0.\n') 
      obj_fichier.write(' 0., 0.1\n') 
      obj_fichier.write('*Surface Interaction, name=IntProp-2\n') 
      obj_fichier.write('1.,\n') 
      obj_fichier.write('*Surface Behavior, pressure-overclosure=HARD\n') 
      obj_fichier.write('**\n') 
      
      # initialisation du champ de température sur tout le domaine
      obj_fichier.write('** PREDEFINED FIELDS\n') 
      obj_fichier.write('**\n') 
      obj_fichier.write('** Name: Predefined Field-1   Type: Temperature\n') 
      obj_fichier.write('*Initial Conditions, type=TEMPERATURE\n') 
      obj_fichier.write('Set-3, 293.\n') 
      obj_fichier.write('**\n') 
      
      # définition des intéractions entres surfaces
      if simulation == 'thermo-meca' :
        obj_fichier.write('** INTERACTIONS\n') 
        obj_fichier.write('**\n') 
        for j in range (0,len(table_Surfaces_interaction),2):
            obj_fichier.write('** Interaction: Int-'+str(int(j/2)+1)+'\n') 
            obj_fichier.write('*Contact Pair, interaction=IntProp-1, type=SURFACE TO SURFACE\n') 
            obj_fichier.write(table_Surfaces_interaction[j][0]+', '+table_Surfaces_interaction[j+1][0]+'\n') 
        obj_fichier.write('**\n')
      
    if ind == 'null' :
       obj_fichier.write(i)

    ligne=i.replace('\n','')
    ligne1=ligne.replace(' ','')
    elements = ligne1.split(',')
    
    # recherche de noms de modules
    if elements[0] =='**PARTS':
      module = 'PARTS'
    elif elements[0] =='*EndPart':
      module = 'null'
    elif elements[0] =='**ASSEMBLY':
      module = 'ASSEMBLY'
    elif elements[0] =='*EndAssembly':
      module = 'null'
    
    # recherche de mots clefs pour modification de la/des lignes suivantes
    if elements[0] =='*Node':
      ind = 'Node'
    elif elements[0] =='*Element':
      ind = 'Element'
    elif elements[0] == '*Nset':
      ind = 'Nset'
    elif elements[0] == '*Elset':
      ind = 'Elset'
    elif elements[0] =='*SolidSection':
      ind = 'Solid Section'
    elif elements[0] =='*SpecificHeat':
      ind = 'Material2'
    elif elements[0] =='*EndInstance':
      ind = 'Set'
    else :
      ind = 'null'

obj_fichier.close()
print(" - Construction du nouveau JDD -> ok")

print("\n-----------------------------------------")
print("                   Fin")
print("-----------------------------------------")


