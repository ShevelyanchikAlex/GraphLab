class Main {
    public static void main(String[] args) {
        WeightedGraph weightedGraph = new WeightedGraph(5);
        weightedGraph.addEdge(0, 1, 1);
        weightedGraph.addEdge(1, 2, 2);
        weightedGraph.addEdge(2, 3, 2);
        weightedGraph.addEdge(2, 4, 4);
        weightedGraph.addEdge(3, 2, 3);
        weightedGraph.addEdge(3, 1, 1);
        weightedGraph.addEdge(4, 3, 5);

        weightedGraph.printAdjMatrix();
        System.out.println();
        weightedGraph.printAdjList();
        System.out.println();
        weightedGraph.searchShortestPath(1, 3);
        weightedGraph.searchShortestPath(2, 3);
        weightedGraph.searchShortestPath(1, 4);
        weightedGraph.searchGraphCenter();

        weightedGraph.findAndPrintAllPaths(1, 3);
    }
}