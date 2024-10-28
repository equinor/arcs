from arcs.setup_functions import GraphGenerator

if __name__ == "__main__":
    g = GraphGenerator(applied_reactions="SCAN_reactions.p")
    g.generatemultidigraph()
    g.save(filename="SCAN_graph.p")
