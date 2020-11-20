import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

class WeightedGraph {

    private int numOfVertex;
    private int[][] adjMatrix;
    private LinkedList<Edge>[] adjList;
    private List<Edge> edgesList;
    private int sizeOfEdgesList;
    private ArrayList<String> listAllPathsStr = new ArrayList<>();


    WeightedGraph(int numOfVertex) {
        this.numOfVertex = numOfVertex;
        this.adjMatrix = new int[numOfVertex][numOfVertex];
        this.adjList = new LinkedList[numOfVertex];
        initAdjList();
        this.edgesList = new ArrayList<Edge>();
        this.sizeOfEdgesList = 0;
    }

    private void initAdjList() {
        for (int i = 0; i < numOfVertex; i++) {
            adjList[i] = new LinkedList<>();
        }
    }


    void addEdge(int source, int destinationVertex, int weight) {
        edgesList.add(new Edge(source, destinationVertex, weight));
        Edge currentEdge = edgesList.get(sizeOfEdgesList);
        addEdgeInAdjMatrixElement(currentEdge);
        addEdgeInAdjList(currentEdge);
        sizeOfEdgesList++;
    }


    private void addEdgeInAdjMatrixElement(Edge currentEdge) {
        adjMatrix[currentEdge.getSourceVertex()][currentEdge.getDestinationVertex()] = currentEdge.getWeight();
    }


    private void addEdgeInAdjList(Edge currentEdge) {
        adjList[currentEdge.getSourceVertex()].addLast(currentEdge);
    }


    void printAdjMatrix() {
        for (int i = 0; i < this.numOfVertex; i++) {
            System.out.print((i + 1) + ":\t");
            for (int j = 0; j < this.numOfVertex; j++) {
                System.out.print(adjMatrix[i][j] + "\t");
            }
            System.out.println();
        }
    }

    void printAdjList() {
        for (int i = 0; i < numOfVertex; i++) {
            LinkedList<Edge> list = adjList[i];
            System.out.print("Source vertex : " + (i + 1));
            for (Edge edge : list) {
                System.out.print(" -> " + (edge.getDestinationVertex() + 1) + " (weight: " + edge.getWeight() + ")");
            }
            System.out.println();
        }
    }


    void searchShortestPath(int sourceVertex, int destinationVertex) {

        double[][] dist = initDistMatrix();
        int[][] next = initNextMatrix();

//        printNextMatrix(next);

        for (int k = 0; k < numOfVertex; k++) {
            for (int i = 0; i < numOfVertex; i++) {
                for (int j = 0; j < numOfVertex; j++) {
                    if (dist[i][k] + dist[k][j] < dist[i][j]) {
                        dist[i][j] = dist[i][k] + dist[k][j];
                        next[i][j] = next[i][k];
                    }
                }
            }
        }

//        printNextMatrix(next);
//        printDistMatrix(dist);
        printPath(dist, next, sourceVertex, destinationVertex);
    }


    private double[][] initDistMatrix() {
        double[][] dist = new double[numOfVertex][numOfVertex];

        for (double[] row : dist)
            Arrays.fill(row, Double.POSITIVE_INFINITY);

        for (int i = 0; i < numOfVertex; i++) {
            for (int j = 0; j < numOfVertex; j++) {
                if (adjMatrix[i][j] != 0) {
                    dist[i][j] = adjMatrix[i][j];
                }
            }
        }
        return dist;
    }

    void searchGraphCenter() {
        double[][] dist = initDistMatrix();
        for (int k = 0; k < numOfVertex; k++) {
            for (int i = 0; i < numOfVertex; i++) {
                for (int j = 0; j < numOfVertex; j++) {
                    if (dist[i][k] + dist[k][j] < dist[i][j]) {
                        dist[i][j] = dist[i][k] + dist[k][j];
                    }
                }
            }
        }

//        printDistMatrix(dist);
        double[] eccentricity = initEccentricityArray(dist);
//        printDistMatrix(dist);
        for (int i = 0; i < numOfVertex; i++) {
            for (int j = 0; j < numOfVertex; j++) {
                if (eccentricity[0] == dist[i][j]) {
                    System.out.println("\nCenter of graph : " + (j + 1));
                    return;
                }
            }
        }
    }


    private double[] initEccentricityArray(double[][] dist) {
        double[] eccentricity = new double[numOfVertex];
        for (int i = 0; i < numOfVertex; i++) {
            if (dist[i][0] != Double.POSITIVE_INFINITY) {
                eccentricity[i] = dist[i][0];
            }
            for (int j = 0; j < numOfVertex; j++) {
                if (dist[i][j] > eccentricity[i] && dist[i][j] != Double.POSITIVE_INFINITY) {
                    eccentricity[i] = dist[i][j];
                }
            }
        }
        Arrays.sort(eccentricity);
        return eccentricity;
    }


    private void printDistMatrix(double[][] dist) {
        System.out.println("Dist matrix :");
        for (double[] doubles : dist) {
            for (int j = 0; j < dist.length; j++) {

                System.out.print(doubles[j] + " ");
            }
            System.out.println();
        }
        System.out.println();
    }


    private int[][] initNextMatrix() {
        int[][] next = new int[numOfVertex][numOfVertex];
        for (int i = 0; i < next.length; i++) {
            for (int j = 0; j < next.length; j++) {
                if (i != j)
                    next[i][j] = j + 1;
            }
        }
        return next;
    }


    private void printNextMatrix(int[][] next) {
        System.out.println("Next matrix :");
        for (int[] ints : next) {
            for (int j = 0; j < next.length; j++) {
                System.out.print(ints[j] + " ");
            }
            System.out.println();
        }
        System.out.println();
    }


     private static void printPath(double[][] dist, int[][] next, int sourceVertex, int destinationVertex) {
        if (sourceVertex != destinationVertex) {
            System.out.println("Pair of vertex\t\tWeight\t\tPath");
            int srcV = sourceVertex;
            int destV = destinationVertex;
            System.out.print(srcV + "->" + destV + "\t\t\t\t" + (int) dist[--sourceVertex][--destinationVertex] + "\t\t\t" + srcV);
            do {
                srcV = next[srcV - 1][destV - 1];
                System.out.print("->" + srcV);
            } while (srcV != destV);
            System.out.println();
        }
    }


    void findAndPrintAllPaths(int sourceVertex, int destinationVertex) {
        boolean[] isVisitedVertex = new boolean[numOfVertex];
        ArrayList<Integer> listAllPaths = new ArrayList<>();

        listAllPaths.add(sourceVertex);
        findAndPrintAllPathsRecursive(sourceVertex, destinationVertex, isVisitedVertex, listAllPaths);
        printAllPaths(sourceVertex, destinationVertex);
    }


    private void findAndPrintAllPathsRecursive(Integer sourceVertex, Integer destinationVertex, boolean[] isVisited, List<Integer> listAllPaths) {

        if (sourceVertex.equals(destinationVertex)) {
            listAllPathsStr.add(String.valueOf(listAllPaths));
            return;
        }

        isVisited[sourceVertex] = true;

        for (Edge i : adjList[sourceVertex]) {
            if (!isVisited[i.getDestinationVertex()]) {

                listAllPaths.add(i.getDestinationVertex());
                findAndPrintAllPathsRecursive(i.getDestinationVertex(), destinationVertex, isVisited, listAllPaths);
                listAllPaths.remove(listAllPaths.size() - 1);
            }
        }
        isVisited[sourceVertex] = false;
    }


    private void printAllPaths(int sourceVertex, int destinationVertex) {
        listAllPathsStr.sort((s1, s2) -> s1.length() - s2.length());
        System.out.println("\n(Here vertex 1 = vertex 0)");
        System.out.println("All paths " + sourceVertex + " -> " + destinationVertex + " :");
        for (String path : listAllPathsStr) {
            printOnePath(path);
        }
    }


    private void printOnePath(String path) {
        for (int i = 0; i < path.length(); i++) {
            if (path.charAt(i) != '[' && path.charAt(i) != ']') {
                if (path.charAt(i) == ',') {
                    System.out.print(" -> ");
                } else {
                    System.out.print(path.charAt(i));
                }
            }
        }
        System.out.println();
    }


    public void setEdgesList(List<Edge> edgesList) {
        this.edgesList = edgesList;
    }

    public List<Edge> getEdgesList() {
        return edgesList;
    }

    public void setSizeOfEdgesList(int sizeOfEdgesList) {
        this.sizeOfEdgesList = sizeOfEdgesList;
    }

    public int getSizeOfEdgesList() {
        return sizeOfEdgesList;
    }

    public void setAdjList(LinkedList<Edge>[] adjList) {
        this.adjList = adjList;
    }

    public LinkedList<Edge>[] getAdjList() {
        return adjList;
    }

    public void setAdjMatrix(int[][] adjMatrix) {
        this.adjMatrix = adjMatrix;
    }

    public int[][] getAdjMatrix() {
        return adjMatrix;
    }

    public int getNumOfVertex() {
        return numOfVertex;
    }

    public void setNumOfVertex(int numOfVertex) {
        this.numOfVertex = numOfVertex;
    }


    static class Edge {
        int sourceVertex;
        int destinationVertex;
        int weight;

        Edge(int sourceVertex, int destinationVertex, int weight) {
            this.sourceVertex = sourceVertex;
            this.destinationVertex = destinationVertex;
            this.weight = weight;
        }

        void setDestinationVertex(int destinationVertex) {
            this.destinationVertex = destinationVertex;
        }

        int getDestinationVertex() {
            return destinationVertex;
        }

        void setSourceVertex(int sourceVertex) {
            this.sourceVertex = sourceVertex;
        }

        int getSourceVertex() {
            return sourceVertex;
        }

        void setWeight(int weight) {
            this.weight = weight;
        }

        int getWeight() {
            return weight;
        }

    }

}
