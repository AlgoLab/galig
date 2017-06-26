from graphviz import Digraph

class Node:
    def __init__(self, label, weight, comp3, comp5):
        self.label = label
        self.weight = weight
        self.comp3 = comp3
        self.comp5 = comp5

    def increment(self, w):
        self.weight+=w

    def decrement(self, w):
        self.weight-=w
        if self.weight < 0:
            self.weight = 0

class SplicingGraph:
    #1-based
    def __init__(self):
        self.labels = []
        self.nodes = {}
        self.edges = {}
        self.new_edges = {}
        self.fake_edges = {}

    def addNode(self, label, w=0, comp3=0, comp5=0):
        if label != "" and label not in self.labels:
            self.labels.append(label)
            self.nodes.update({len(self.labels):Node(label, w, comp3, comp5)})

    def addEdge(self, n1, n2, t, w=0):
        if t == 'e':
            if (n1,n2) not in self.edges:
                self.edges.update({(n1,n2):w})
        elif t == 'n':
            if (n1,n2) not in self.new_edges:
                self.new_edges.update({(n1,n2):w})
        elif t == 'f':
            if (n1,n2) not in self.fake_edges:
                self.fake_edges.update({(n1,n2):w})

    def incrementNode(self, n1, w=1):
        if n1 in self.nodes:
            self.nodes[n1].increment(w)
            #print("Nodes {} incremented".format(self.nodes[n1-1]))

    def decrementNode(self, n1, w=1):
        if n1 in self.nodes:
            self.nodes[n1].decrement(w)
            #print("Nodes {} decremented".format(self.nodes[n1-1]))

    def incrementEdge(self, n1, n2, w=1):
        if (n1,n2) in self.edges:
            self.edges[(n1,n2)] += w
        elif (n1,n2) in self.new_edges:
            self.new_edges[(n1,n2)] += w
        elif (n1,n2) in self.fake_edges:
            self.fake_edges[(n1,n2)] += w

    def decrementEdge(self, n1, n2, w=1):
        if (n1,n2) in self.edges:
            if self.edges[(n1,n2)] > 0:
                self.edges[(n1,n2)] -= w
        elif (n1,n2) in self.new_edges:
            if self.new_edges[(n1,n2)] > 0:
                self.new_edges[(n1,n2)] -= w
        elif (n1,n2) in self.fake_edges:
            if self.fake_edges[(n1,n2)] > 0:
                self.fake_edges[(n1,n2)] -= w

    def getLabel(self, i):
        return self.nodes[i].label

    def getIndex(self, label):
        #Not handled exception!!!
        return self.labels.index(label)+1

    def clean(self, nodes_min, edges_min):
        edges_to_remove = []
        for label in self.labels:
            if label == "":
                continue
            index = self.labels.index(label)+1
            has_fake = False
            if self.nodes[index].weight <= nodes_min:
                for (n1,n2),cov in self.fake_edges.items():
                    if index == n1 or index == n2:
                        has_fake = True
                if has_fake:
                    continue
                self.labels[index-1] = ""
                self.nodes.pop(index)
                for (n1,n2),cov in self.edges.items():
                    if index == n1 or index == n2:
                        edges_to_remove.append((n1,n2))
                for (n1,n2),cov in self.new_edges.items():
                    if index == n1 or index == n2:
                        edges_to_remove.append((n1,n2))

        for (n1,n2),cov in self.edges.items():
            if cov <= edges_min:
                edges_to_remove.append((n1,n2))
        for (n1,n2),cov in self.new_edges.items():
            if cov <= edges_min:
                edges_to_remove.append((n1,n2))
        for edge in edges_to_remove:
            try:
                self.edges.pop(edge)
                continue
            except KeyError:
                pass
            try:
                self.new_edges.pop(edge)
                continue
            except KeyError:
                pass

    def print(self):
        print("### NODES ###")
        for label in self.labels:
            if label != "":
                index = self.labels.index(label)+1
                print("{} {}: {}".format(index, label, self.nodes[index].weight))
        print("### EDGES ###")
        for (n1,n2),c in self.edges.items():
            print("{}->{}: {}".format(n1, n2, c))
        print("### NEW EDGES ###")
        for (n1,n2),c in self.new_edges.items():
           print("{}->{}: {}".format(n1, n2, c))
        print("### FAKE EDGES ###")
        for (n1,n2),c in self.fake_edges.items():
           print("{}->{}: {}".format(n1, n2, c))

    def save(self, name):
        g = Digraph('G', filename="./{}.gv".format(name))#os.path.join(OUT, "graph.gv"))
        g.attr('node', shape='circle')
        for label in self.labels:
            if label != "":
                index = self.labels.index(label)
                n_label = "{} ({})".format(label, self.nodes[index+1].weight)
                g.node(n_label)
        for (n1,n2),cov in self.edges.items():
            n1_label = "{} ({})".format(self.labels[n1-1], self.nodes[n1].weight)
            n2_label = "{} ({})".format(self.labels[n2-1], self.nodes[n2].weight)
            g.edge(n1_label, n2_label, label=str(cov), color = "black")
        for (n1,n2),cov in self.new_edges.items():
            n1_label = "{} ({})".format(self.labels[n1-1], self.nodes[n1].weight)
            n2_label = "{} ({})".format(self.labels[n2-1], self.nodes[n2].weight)
            g.edge(n1_label, n2_label, label=str(cov), color = "red")
        for (n1,n2),cov in self.fake_edges.items():
            n1_label = "{} ({})".format(self.labels[n1-1], self.nodes[n1].weight)
            n2_label = "{} ({})".format(self.labels[n2-1], self.nodes[n2].weight)
            g.edge(n1_label, n2_label, label=str(cov), color = "black", style="dashed")
        g.render()
