
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 23 15:00:53 2018

@author: haswani
"""
import math
from skan import csr
import numpy as np
from skimage.morphology import skeletonize,skeletonize_3d
from astropy.io import fits
from skimage.morphology import binary_opening,binary_dilation,binary_erosion,closing,dilation,opening,ball,remove_small_holes
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import networkx as nx

# import imageio as iio
# from skan import draw
#import networkx as nx

def returnSubgraph(G):
    '''
    Return list of all the subgraph of networkX graph G
    
    Parameters
    ----------
    G : NetworkX Graph 

   
    
    
    Returns
    -------
    List containing node information of all the subgraphs of the Graph G
    '''
    #requires Networks as nx
    sub_graphs = nx.connected_component_subgraphs(G)

    subgraphlist = []

    for i, sg in enumerate(sub_graphs):
    #print("subgraph {} has {} nodes".format(i, sg.number_of_nodes()))
    #print("\tNodes:", sg.nodes())
        subgraphlist.append(sg.nodes())
    
    return subgraphlist

def returnEndpointsAndJunction(H3):
    '''
    Return list of all the End Points and Junction nodes of networkX graph
    
    Parameters
    ----------
    H3 : NetworkX Graph 

   
    
    
    Returns
    -------
    List containing node information of all the endpoints and junction nodes
    of the Graph G
    '''
    #requires Networks as nx
    endPoints = []
    junction = []
    for n in H3.nodes:
        if H3.degree(n)==1:
            endPoints.append(n)
        elif H3.degree(n)>2:
            junction.append(n)
    return endPoints,junction
            

def network_plot_3D(G, angle, save=False):
    '''
    Gives a 3D plot for networkX using coordinates information of the nodes
    
    Parameters
    ----------
    G : NetworkX Graph
    angle : Angel to view the graph plot
    save : boolen value when true saves the plot

   
    
    
    Returns
    -------
    Plot for networkX graph.
    '''
    # Get node positions
    pos = nx.get_node_attributes(G, 'pos')

    # Get number of nodes

    # Get the maximum number of edges adjacent to a single node
    edge_max = max([G.degree(i) for i in G.nodes])

    # Define color range proportional to number of edges adjacent to a single node
    colors = [plt.cm.plasma(G.degree(i)/edge_max) for i in G.nodes]

    # 3D network plot
    with plt.style.context(('ggplot')):

        fig = plt.figure(figsize=(10,7))
        ax = Axes3D(fig)

        # Loop on the pos dictionary to extract the x,y,z coordinates of each node
        for i, (key, value) in enumerate(pos.items()):
            xi = value[0]
            yi = value[1]
            zi = value[2]

            # Scatter plot
            ax.scatter(xi, yi, zi, c=colors[i], s=20+20*G.degree(key), edgecolors='k', alpha=0.7)

        # Loop on the list of edges to get the x,y,z, coordinates of the connected nodes
        # Those two points are the extrema of the line to be plotted
        for i,j in enumerate(G.edges()):

            x = np.array((pos[j[0]][0], pos[j[1]][0]))
            y = np.array((pos[j[0]][1], pos[j[1]][1]))
            z = np.array((pos[j[0]][2], pos[j[1]][2]))

        # Plot the connecting lines
            ax.plot(x, y, z, c='black', alpha=0.5)

    # Set the initial view
    ax.view_init(30, angle)

    # Hide the axes
    #ax.set_axis_off()

    if save is not False:
        plt.savefig(str(angle).zfill(3)+".png")
        plt.close('all')
    else:
        plt.show()

    return



#cube = fits.getdata('ngc4321_co21_12m+7m+tp_pbcorr_round_k_correct_mask.fits')
##cube = fits.getdata('ngc3627_co21_12m+7m+tp_mask.fits')
#parameter = 5


def skeletonize3D(cube):
    """
    Return the 3D skeleton dkel2 and converted it in form of graph and return
    pixel graph as a SciPy CSR matrix in which entry (i,j) is 0 if pixels 
    i and j are not connected, and otherwise is equal to the distance between
    pixels i and j in the skeleton.

    coordinates (in pixel units) of the points in the pixel graph. Finally, degrees 
    is an image of the skeleton, with each skeleton pixel containing the number of neighboring pixels.

 

    Parameters
    ----------
    cube : 3d numpy array

   
    
    
    Returns
    -------
    dskel2 :3D skeletonized binary image(numpy array)
    Pixel Graph,
    coordinates,
    degree 
       
    """
    
    selem = ball(3)

    ddilate = dilation(cube, selem)
    dclose = closing(ddilate)
    dskel2 = skeletonize_3d(dclose)


    pixel_graph0, coordinates0, degrees0 = csr.skeleton_to_csgraph(dskel2)

    return dskel2, pixel_graph0, coordinates0, degrees0    

    

def pruning3D(pixel_graph0,coordinates0,parameter):
    '''
    
    Convert pixel_Graph into Networkx graph and prune all the branches and 
    subgraph less then parameter length
    
    Parameters :
         pixel_graph0: is a SciPy CSR matrix in which entry (i,j) is 0 if pixels 
    i and j are not connected, and otherwise is equal to the distance between
    pixels i and j in the skeleton.

    coordinates0 : Coordinates (in pixel units) of the points in the pixel graph.
           
    parameter : pruning all the branches shorter then this length and 
                all subgraph less then this length are removed      
    
    '''        
        
    cutnodes = []
    G = nx.from_scipy_sparse_matrix(pixel_graph0)

    for node in G.node:
        G.node[node]['pos'] = coordinates0[node]

    print("Number of nodes before Pruning:")
    print(nx.number_of_nodes(G))


# extract subgraphs
    sub_graphs = nx.connected_component_subgraphs(G)

    subgraphlist = []

    for i, sg in enumerate(sub_graphs):
    #print("subgraph {} has {} nodes".format(i, sg.number_of_nodes()))
    #print("\tNodes:", sg.nodes())
        subgraphlist.append(sg.nodes())
        if i == 3: # change this value to check with different subgraphs
            L = sg.nodes()
    #print("\tEdges:", sg.edges())

#H = nx.Graph(G.subgraph(L))  - this is working

    for j in subgraphlist:
        H = nx.Graph(G.subgraph(j)) # for all the subgraphs
#print(list(itertools.chain(H.nx.Graph.edges_iter())))

#for n in H.nodes:
    #    print(G.degree(n))
        endPoints = []

        for n in H.nodes:
            if H.degree(n)==1:
                endPoints.append(n)
    #print(endPoints) #all the endpoints

        junction = []

        for n in H.nodes:
            if H.degree(n)>2:
                junction.append(n) # all the junction points


        for k in junction: #for all the junction nodes removing the branches
#            p = nx.shortest_path(H,source=k,target=endPoints[0])
#            print("Path :")
#            print(len(p))
#
          #  minimum = 30 #just a large number
                for i in endPoints:
                    if nx.has_path(H,k,i) and H.degree(k) > 2:

                        p = nx.shortest_path(H,source=k,target=i)
          #          print("Path :")
           #         print(i)
           #         print(len(p))
                        length = len(p)
                        if length < parameter:
                        #minimum = length
                            cutnode = p[1]
            #            print("TO cut: ")
             #           print(cutnode)


                            if H.has_edge(k, cutnode):
                                H.remove_edge(k,cutnode)
              #              print("SubGraph Number")
               #             print(j)
                #            print("Cutting edges ")
                 #           print(k, cutnode)
                                cutnodes.append(cutnode)
     #Removing Subgraph less then parameter length(here 5nodes)

    
                for n in cutnodes:
                    if G.has_edge(k,n):
                        G.remove_edge(k,n)
                        print("removed : ",k,"---",n)
            #print(cutnodes)



#Removing nodes
        if nx.number_of_nodes(H) <= parameter:
            print("Removing nodes : ",H.nodes())
            print(H.nodes())
            G.remove_nodes_from(H.nodes())

    print("Number of nodes After removing Subgraphs < parameter length:")
    print(nx.number_of_nodes(G))

#removing Subgraphs afer removing edges

    sub_graphs2 = nx.connected_component_subgraphs(G)

    subgraphlist2 = []

    for i, sg in enumerate(sub_graphs2):
    #print("subgraph {} has {} nodes".format(i, sg.number_of_nodes()))
    #print("\tNodes:", sg.nodes())
        subgraphlist2.append(sg.nodes())

    for j in subgraphlist2:
        H2 = nx.Graph(G.subgraph(j))
        if nx.number_of_nodes(H2) <= parameter:
            print("removing nodes 2")
            print(H.nodes())
            G.remove_nodes_from(H2.nodes())


    print("Number of nodes After Pruning:")
    print(nx.number_of_nodes(G))

    return G

######################################################
############TEST DATA 1 without Junction###############################
def testGraph1():
    '''
    return Graph to check Test case when there are no junction
    '''
    Gtest=nx.Graph()

    Gtest.add_node(1, pos=[0., 0., 0.])
    Gtest.add_node(2, pos=[0., 0., 1.])
    Gtest.add_node(3, pos=[0., 0., 2.])
    Gtest.add_node(4, pos=[0., 0., 3.])
    Gtest.add_node(5, pos=[1., 0., 3.])
    Gtest.add_node(6, pos=[2., 0., 3.])
    Gtest.add_node(7, pos=[3., 0., 3.])
    Gtest.add_node(8, pos=[4., 0., 3.])

    Gtest.add_edges_from([[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,8]])
    
    for node in Gtest.node:
        print(Gtest.node[node]['pos'])
    return Gtest    


#####################################################


######################################################
############TEST DATA 2 with Junction###############################
    

def testGraph2():
    '''
    return Graph to check Test case when there are junction
    '''
    Gtest2=nx.Graph()

    Gtest2.add_node(1, pos=[0., 0., 0.])
    Gtest2.add_node(2, pos=[0., 0., 1.])
    Gtest2.add_node(3, pos=[0., 0., 2.])
    Gtest2.add_node(4, pos=[0., 0., 3.])
    Gtest2.add_node(5, pos=[1., 0., 3.])
    Gtest2.add_node(6, pos=[2., 0., 3.])
    Gtest2.add_node(7, pos=[0., 0., 3.])
    Gtest2.add_node(8, pos=[0., 0., 4.])
    
    Gtest2.add_edges_from([[1,2],[2,3],[3,4],[4,5],[5,6],[4,7],[7,8]])

    for node in Gtest2.node:
        print(Gtest2.node[node]['pos'])
    return Gtest2    


def plotAndReturnLongestFilaments(G):
    finalGraph = nx.Graph()
    node_number = []
    cordinates = []
    spatial_length = []
    veloctiy_length = []
#will zip them all together in end

#sub_graphs3 = nx.connected_component_subgraphs(G)
#checking test case
    sub_graphs3 = nx.connected_component_subgraphs(G)

#extracting all subgraphs to work on each subgraphs individually

    subgraphlist3 = []

    for i, sg in enumerate(sub_graphs3):
    #print("subgraph {} has {} nodes".format(i, sg.number_of_nodes()))
    #print("\tNodes:", sg.nodes())
        subgraphlist3.append(sg.nodes())
     #   if i == 3: # change this value to check with different subgraphs
  #          L = sg.nodes()

#subgraphlist3 is the required list of all subgraphs

    resultGraph = nx.Graph()
    longestFilamentGraph = nx.Graph()
    for j in subgraphlist3:
        H3 = nx.Graph(G.subgraph(j))
        endPoints = [] #endpoints of each subgraph
        junction = []   #junction of each subgraph
        for n in H3.nodes:
            if H3.degree(n)==1:
                endPoints.append(n)
            elif H3.degree(n)>2:
                junction.append(n)
        
        if len(endPoints) == 0:
            raise ValueError("Found no end points.")

        #print("EndPoints: ")
        #print(endPoints)

        used_ends = []

        for end in endPoints:
         
        # Check whether any branches need to be traversed
            if len(list(set(endPoints) - set(used_ends))) == 0:
                break

            used_ends.append(end)

        # Iterate through nodes until hitting a junction or end point
            spat_length = 0.
            vel_length = 0.

            curr_node = end
            used_nodes = []
        
        ##
        #    print("EndPoint coordinate")
        
        #cordinates of the endpoint
            pointer_node = H3.node[curr_node]['pos']
        #    print(pointer_node)
                

            
        # Traverse along the path until hitting a junction or end
            while True:
         #       print("In loop")
                neighbors = list(H3.neighbors(curr_node))
            
            # Find the new node
                if len(neighbors) > 1:
                # Remove previously visited nodes
                    neighbors = list(set(neighbors) - set(used_nodes))

   
                neighbour = neighbors[0]
                next_node = H3.node[neighbour]['pos']
           #     print("Neighbour Pixel ")
            
            # Check for stopping criterion if the next node is an end point
            # or junction
                if neighbour in endPoints:
                # Setup to skip that end point
                
                #Make it a function Laterr
                #calculate the spatial length
                    res1 = abs(next_node[1]-pointer_node[1])
                    res2 = abs(next_node[2]-pointer_node[2])
                    x = math.pow(res1, 2) + math.pow(res2, 2)
                    y = math.sqrt(x)
                    spat_length = spat_length + y 
                
                #calculating velocity length 
                    vel_length = vel_length + abs(next_node[0]-pointer_node[0])
                
                
                    used_ends.append(neighbour)
                
                    resultGraph.add_node(neighbour) #graph to find the longest Filament
                    resultGraph.node[neighbour]['pos'] = H3.node[neighbour]['pos']
                    resultGraph.add_edge(end,neighbour,S_lenght=spat_length, V_lenght = vel_length)
                
                ##
                    finalGraph.add_node(neighbour) #graph to find the longest Filament
                    finalGraph.node[neighbour]['pos'] = H3.node[neighbour]['pos']
                    finalGraph.add_edge(end,neighbour,S_lenght=spat_length, V_lenght = vel_length)
                ##
                    break
                elif neighbour in junction:
                # End on junctions
    
                #calculate the spatial length
                    res1 = abs(next_node[1]-pointer_node[1])
                    res2 = abs(next_node[2]-pointer_node[2])
                    x = math.pow(res1, 2) + math.pow(res2, 2)
                    y = math.sqrt(x)
                    spat_length = spat_length + y 
                
                #calculating velocity length 
                    vel_length = vel_length + abs(next_node[0]-pointer_node[0])
                
                
                    resultGraph.add_node(neighbour) #graph to find the longest Filament
                    resultGraph.node[neighbour]['pos'] = H3.node[neighbour]['pos']
                    resultGraph.add_edge(end,neighbour,S_length=spat_length, V_length = vel_length)
                
                ##
                    finalGraph.add_node(neighbour) #graph to find the longest Filament
                    finalGraph.node[neighbour]['pos'] = H3.node[neighbour]['pos']
                    finalGraph.add_edge(end,neighbour,S_length=spat_length, V_length = vel_length)
                ##
                    break
                else:
                # Update along a branch
                    used_nodes.append(curr_node)
                
                #calculate the spatial length
                    res1 = abs(next_node[1]-pointer_node[1])
                    res2 = abs(next_node[2]-pointer_node[2])
                    x = math.pow(res1, 2) + math.pow(res2, 2)
                    y = math.sqrt(x)
                    spat_length = spat_length + y 
                
                #calculating velocity length 
                    vel_length = vel_length + abs(next_node[0]-pointer_node[0])
                
                    curr_node = neighbour
                #pointer Node contain cordinates of the updated current_node which is the neighbour node
                    pointer_node = H3.node[curr_node]['pos']
     #               print(pointer_node)
    #                print(vel_length)
   #                 print(spat_length)
        #updating List
            node_number.append(end)
            resultGraph.add_node(end) #graph to find the longest Filament
        
            finalGraph.add_node(end) ##
        
        
        
            cordinates.append(H3.node[end]['pos'])
        
            finalGraph.node[end]['pos'] = H3.node[end]['pos'] ##
        
            resultGraph.node[end]['pos'] = H3.node[end]['pos'] #adding cordinate attribute to Result Graph
        
            spatial_length.append(spat_length)
            veloctiy_length.append(vel_length)        
        
        longestLength = 0 

    
        e_p,j_p = returnEndpointsAndJunction(finalGraph)
        for end1 in e_p:
            for end2 in e_p:
                if nx.has_path(finalGraph,end1,end2):
                    path = nx.shortest_path(finalGraph, source=end1, target=end2, weight='S_length')
                    pathLength= nx.shortest_path_length(finalGraph,source=end1, target = end2 , weight = 'S_length')
                    if pathLength > longestLength:
                        longestLength = pathLength
                        requiredPath = path
                        start = end1 # will be used in the original graph to get the longest filaments
                        finish = end2
                
            #here can also return endPoints
        #    print("Longest Path with Spatial Length as edge attribute",longestLength)
     #   print("Requierd Path", requiredPath)       
        longestFilament = nx.shortest_path(H3,source=start,target=finish)
        print("Longest Path for Filament" , longestFilament)
        finalGraph.clear()
        temp_graph = H3.subgraph(longestFilament)
    #nx.draw(temp_graph,label=True)
    #network_plot_3D(temp_graph,30, save=False)
        longestFilamentGraph = nx.union(longestFilamentGraph,temp_graph)
    #uncomment below code to print more details
    
    #nx.draw(longestFilamentGraph,label=True)
    #network_plot_3D(longestFilamentGraph,30, save=False)
    #result = list(zip(node_number, cordinates, spatial_length, veloctiy_length))


# Converting itertor to set
    #print(result)

    #print("to Check:")

    #nx.draw(finalGraph, with_labels=True)
    #for node in resultGraph.node:
    #    print(resultGraph.node[node]['pos'])
    #
    #print(resultGraph.edges(data=True))
    #for node in longestFilamentGraph.node:
    #   print(longestFilamentGraph.node[node]['pos'])

    #print(longestFilamentGraph.edges(data=True))


    network_plot_3D(longestFilamentGraph,30, save=False)
    return longestFilamentGraph

    
    ##############################
def filVelocityAndSpatialLength3D(Gtest2):
    '''
    Gtest2 Networkx graph representation of filaments is used to find 
    spatial length and velocity length of the filaments
    
    Parameters
    ----------
    Gtest2 : 3D Networkx Graph 

   
    
    
    Returns
    -------
    List containing : node number of end Points, Coordinates of this end point,
    Spatial Length (of this end point to the next junction or endpoint),
    Velocity length (of this end point to the next junction or endpoint)
    
    '''


    node_number = []
    cordinates = []
    spatial_length = []
    veloctiy_length = []
    #will zip them all together in end
    
    #print(argh)
    
    #sub_graphs3 = nx.connected_component_subgraphs(G)
    #checking test case
    sub_graphs3 = nx.connected_component_subgraphs(Gtest2)
    
    #extracting all subgraphs to work on each subgraphs individually
    
    subgraphlist3 = []
    
    for i, sg in enumerate(sub_graphs3):
        #print("subgraph {} has {} nodes".format(i, sg.number_of_nodes()))
        #print("\tNodes:", sg.nodes())
        subgraphlist3.append(sg.nodes())
    #    if i == 3: # change this value to check with different subgraphs
    #        L = sg.nodes()
        
    #subgraphlist3 is the required list of all subgraphs
        
    
    for j in subgraphlist3:
        H3 = nx.Graph(Gtest2.subgraph(j))
        endPoints = [] #endpoints of each subgraph
        junction = []   #junction of each subgraph
        for n in H3.nodes:
            if H3.degree(n)==1:
                endPoints.append(n)
            elif H3.degree(n)>2:
                junction.append(n)
            
        if len(endPoints) == 0:
            raise ValueError("Found no end points.")
        
    
     
        print("EndPoints: ")
        print(endPoints)

        used_ends = []

        for end in endPoints:
            print("End Node Label")
            print(end)
        # Check whether any branches need to be traversed
            if len(list(set(endPoints) - set(used_ends))) == 0:
                break

            used_ends.append(end)

        # Iterate through nodes until hitting a junction or end point
            spat_length = 0.
            vel_length = 0.

            curr_node = end
            used_nodes = []
        
        
            print("EndPoint coordinate")
        
        #cordinates of the endpoint
            pointer_node = H3.node[curr_node]['pos']
            print(pointer_node)
                

            
        # Traverse along the path until hitting a junction or end
            while True:
                
                neighbors = list(H3.neighbors(curr_node))
            
            # Find the new node
                if len(neighbors) > 1:
                # Remove previously visited nodes
                    neighbors = list(set(neighbors) - set(used_nodes))

                neighbour = neighbors[0]
                next_node = H3.node[neighbour]['pos']
            
            

            # stopping criterion if the next node is an end point
            # or junction
                if neighbour in endPoints:
                # Setup to skip that end point
            
                #calculate the spatial length
                    res1 = abs(next_node[1]-pointer_node[1])
                    res2 = abs(next_node[2]-pointer_node[2])
                    x = math.pow(res1, 2) + math.pow(res2, 2)
                    y = math.sqrt(x)
                    spat_length = spat_length + y 
                
                #calculating velocity length 
                    vel_length = vel_length + abs(next_node[0]-pointer_node[0])
                
                
                    used_ends.append(neighbour)
                    
                    break
                elif neighbour in junction:
                # End on junctions
                                #Make it a function Laterr

                #calculate the spatial length
                    res1 = abs(next_node[1]-pointer_node[1])
                    res2 = abs(next_node[2]-pointer_node[2])
                    x = math.pow(res1, 2) + math.pow(res2, 2)
                    y = math.sqrt(x)
                    spat_length = spat_length + y 
                
                #calculating velocity length 
                    vel_length = vel_length + abs(next_node[0]-pointer_node[0])
                
                    print("Found Junction")
                    break
                else:
                # Update along a branch
                    used_nodes.append(curr_node)
                
                #calculate the spatial length
                    res1 = abs(next_node[1]-pointer_node[1])
                    res2 = abs(next_node[2]-pointer_node[2])
                    x = math.pow(res1, 2) + math.pow(res2, 2)
                    y = math.sqrt(x)
                    spat_length = spat_length + y 
                
                #calculating velocity length 
                    vel_length = vel_length + abs(next_node[0]-pointer_node[0])
                
                    curr_node = neighbour
                #pointer Node contain cordinates of the updated current_node which is the neighbour node
                    pointer_node = H3.node[curr_node]['pos']
                    print(pointer_node)
                    print(vel_length)
                    print(spat_length)
        #updating List
            node_number.append(end)
            cordinates.append(H3.node[end]['pos'])
            spatial_length.append(spat_length)
            veloctiy_length.append(vel_length)        
    

    result = list(zip(node_number, cordinates, spatial_length, veloctiy_length))

    return result

def plot3dnumpy(numpy3d):
    '''
    Plots numpy3d using metplotlib
    
    Parameters
    ----------
    numpy3d : 3D numpy array 

    '''
    
    out = np.where(numpy3d)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(out[0],out[1],out[2],'r,')
    return    



