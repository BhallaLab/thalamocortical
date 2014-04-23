import igraph as ig
from trbnetdata import TraubFullNetData
tndata = TraubFullNetData()
g = ig.Graph(0)
g.add_vertices(len(tndata.celltype))
g.vs['label'] = tndata.celltype
g.vs['count'] = tndata.cellcount
g.vs['ectopicinterval'] = tndata.ectopic_interval
g.write_graphml('unconnected_celltype_net.graphml')
