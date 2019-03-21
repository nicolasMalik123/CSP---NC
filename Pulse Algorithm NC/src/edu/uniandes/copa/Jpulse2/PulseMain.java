/**
 * This is the main class for the pulse algorithm. - This pulse includes path completion strategy
 * 
 * Ref.: Lozano, L. and Medaglia, A. L. (2013). 
 * On an exact method for the constrained shortest path problem. Computers & Operations Research. 40 (1):378-384.
 * DOI: http://dx.doi.org/10.1016/j.cor.2012.07.008 
 * 
 * 
 * @author L. Lozano & D. Duque & N. Cabrera
 * @affiliation Universidad de los Andes - Centro para la Optimización y Probabilidad Aplicada (COPA)
 * @url http://copa.uniandes.edu.co/
 * 
 */

package edu.uniandes.copa.Jpulse2;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.StringTokenizer;

public class PulseMain {

	/**
	 * This method runs the pulse
	 * @param args
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static void main(String[] args) throws IOException, InterruptedException {
		
		//Create the test instance:
		
			File testFile = new File("./networks/Config.txt");
			
			BufferedReader bufRedr = new BufferedReader(new FileReader(testFile));
			
			String actLine = null;
			
			String [] information = new String [6];
			
			int rowA = 0;
			int colA = 0;
			
			while((actLine = bufRedr.readLine()) != null && rowA < 6){	
				String [] info = actLine.split(":");
				information[rowA] = info[1];
				rowA++;
			}

		//Para la dirección backward: 
		
	
			DataHandler dataA = new DataHandler(Integer.parseInt(information[2]),Integer.parseInt(information[1]),Integer.parseInt(information[5]),Integer.parseInt(information[4]),1,information[0]);	
			dataA.ReadDimacsB();
			PulseGraph network = createGraph(dataA);	
			network.SetConstraint(Integer.parseInt(information[3]));
			
		//Para la dirección forward: 
		
			/**
			DataHandler dataA = new DataHandler(Integer.parseInt(information[2]),Integer.parseInt(information[1]),Integer.parseInt(information[4]),Integer.parseInt(information[5]),1,information[0]);
			dataA.ReadDimacsF();
			PulseGraph network = createGraph(dataA);	
			network.SetConstraint(Integer.parseInt(information[3]));
			*/
		
		//Initial bounds
		
			SP(dataA,network);
		
		//Pulse algorithm parameters:
		
			//First primal bound
		
			int MD=network.getVertexByID(dataA.Source-1).getMaxDist();
			network.setPrimalBound(MD);
			network.TimeStar = network.getVertexByID(dataA.Source-1).getMinTime();

			//Size of Q
			
			dataA.numLabels = 10;

			//Initial weights
			
			int[] initialWeights = new int[2];
			initialWeights[0] = 0;
			initialWeights[1] = 0;
		
			//Initial path
			
			ArrayList<Integer> finalPath = new ArrayList();
			
		//Run the pulse and recover the path:
		
			//Starts the clock
			
			double iniTime = System.nanoTime(); 
		
			network.getVertexByID(dataA.Source-1).pulse(initialWeights, finalPath);
			finalPath = completeThePath(dataA,network);
			
			//Ends the clock
			double finalTime = (System.nanoTime() - iniTime)/1000000000;
		
		//Print results
		
		System.out.println("The final path is: "+finalPath);
		System.out.println("The total cost is: "+network.PrimalBound);
		System.out.println("The total time is: "+network.TimeStar);
		System.out.println("The time limit is: "+network.TimeC);
		System.out.println("The computational time is: "+finalTime+" s");
	}
	
	/**
	 * This method creates a network
	 * @param data
	 * @return the final graph
	 */
	private static PulseGraph createGraph(DataHandler data) {
		int numNodes = data.NumNodes;
		PulseGraph Gd = new PulseGraph(numNodes);
		for (int i = 0; i < numNodes; i++) {
			if(i!=(data.LastNode-1)){
				Gd.addVertex(new VertexPulse(i) ); //Primero lo creo, y luego lo meto. El id corresponde al número del nodo
			}
		}
		FinalVertexPulse vv = new FinalVertexPulse(data.LastNode-1);
		Gd.addFinalVertex(vv);
		for(int i = 0; i <data.NumArcs; i ++){
			Gd.addWeightedEdge( Gd.getVertexByID(data.Arcs[i][0]), Gd.getVertexByID(data.Arcs[i][1]),data.Distance[i],data.Time[i], i);			
		}
		return Gd;
	}
	
	/**
	 * This method runs a shortest path for cost and times
	 * @param data
	 * @param network
	 * @throws InterruptedException
	 */
	private static void SP(DataHandler data, PulseGraph network) throws InterruptedException {
		// Create two threads and run parallel SP for the initialization		
		Thread tTime = new Thread();
		Thread tDist = new Thread();
		// Reverse the network and run SP for distance and time 
		DukqstraDist spDist = new DukqstraDist(network, data.LastNode-1);
		DukqstraTime spTime = new DukqstraTime(network, data.LastNode-1);
		tDist = new Thread(new ShortestPathTask(1, spDist, null));
		tTime = new Thread(new ShortestPathTask(0, null,  spTime));
		tDist.start();
		tTime.start();
		tDist.join();
		tTime.join();
	}
	
	/**
	 * This method completes the path using the minimum cost or minimum time path
	 * @param data
	 * @param network
	 * @return the final path
	 */
	public static ArrayList<Integer> completeThePath(DataHandler data, PulseGraph network){
		ArrayList<Integer> path = new ArrayList<Integer>();
		for(int i = 0;i<network.Path.size();i++){
			path.add(network.Path.get(i));
		}
		
		
		int nodoInicial = network.getFinalNode();
		boolean termine = false;
		double costoAcumulado = network.getFinalCost();
		double tiempoAcumulado = network.getFinalTime();
		
		
		while(termine == false) {
			int nodoActual = data.LastNode;
			for(int i = 0; i < network.getVertexByID(nodoInicial).magicIndex.size(); i++) {
				int e = (Integer) network.getVertexByID(nodoInicial).magicIndex.get(i);
				int a = DataHandler.Arcs[e][1];
				
				if(network.best == 1){
					if(costoAcumulado + DataHandler.Distance[e] + network.getVertexByID(a).minDist == network.PrimalBound ) {
						if(tiempoAcumulado+ DataHandler.Time[e] + network.getVertexByID(a).maxTime == network.TimeStar) {
							costoAcumulado+=DataHandler.Distance[e];
							tiempoAcumulado+=DataHandler.Time[e];
							nodoActual = a;	
						}
					}
				}else{
					if(costoAcumulado + DataHandler.Distance[e] + network.getVertexByID(a).maxDist == network.PrimalBound ) {
						if(tiempoAcumulado+ DataHandler.Time[e] + network.getVertexByID(a).minTime == network.TimeStar) {
							costoAcumulado+=DataHandler.Distance[e];
							tiempoAcumulado+=DataHandler.Time[e];
							nodoActual = a;	
						}
					}
				}
			}
		
			path.add(nodoActual);
			if(nodoActual == data.LastNode) {
				termine = true;
			}else {
				nodoInicial = nodoActual;
			}
		}
		
		return path;
	}
	
	/**
	 * This method tells if the final path was found by a minimum time or a minimum cost path completion
	 * @param network
	 * @return 1: cost 2: time
	 */
	public static String whoFindThePath(PulseGraph network) {
		String rta = "There is no point";
		if(network.best == 1) {
			return "Forward cost path completion";
		}
		else if (network.best == 2) {
			return "Forward time path completion";
		}
		return rta;
	}

	

	
}
