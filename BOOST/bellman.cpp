#include <boost/config.hpp>
#include <vector>
#include <iostream>
#include <assert.h>
#include <queue>
#include <utility>
#include <boost/random.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/lockfree/queue.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/bellman_ford_shortest_paths.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <LEDA/graph/graph.h>
#include <LEDA/graph/edge_array.h>
#include <LEDA/graph/node_array.h>
#include <LEDA/graph/graph_misc.h>
#include <LEDA/core/random_source.h>
#include <LEDA/graph/templates/shortest_path.h>

#define nullVertex -1

using namespace boost;
typedef leda::node node;
typedef leda::edge edge;
struct NodeInfo
{
    unsigned int distance;
    bool visited;
};
struct EdgeProperties
{
    int weight;
};

// typedef adjacency_list<vecS, vecS, bidirectionalS, no_property, EdgeProperties> Graph;
typedef boost::adjacency_list<boost::setS, boost::vecS, boost::bidirectionalS, NodeInfo, EdgeProperties> Graph;

typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef property_map<Graph, int EdgeProperties::*>::type weight_pmap;
typedef std::pair<int, int> E;

typedef leda::graph LGraph;
using leda::edge_array;
using leda::node_array;
typedef boost::graph_traits<Graph>::vertex_descriptor vertex;
typedef boost::graph_traits<Graph>::edge_descriptor arc;
typedef boost::graph_traits<Graph>::edge_iterator arcIt;
typedef boost::graph_traits<Graph>::vertex_iterator vertexIt;
typedef boost::graph_traits<Graph>::out_edge_iterator outArcIt;
typedef boost::graph_traits<Graph>::adjacency_iterator aIt;

typedef boost::property_map<Graph, int EdgeProperties::*>::type WeightPrMap;
typedef boost::property_map<Graph, unsigned int NodeInfo::*>::type DistanceMap;
typedef boost::property_map<Graph, bool NodeInfo::*>::type VisitedMap;

class MyVisitor : public boost::default_dfs_visitor
{
public:
    MyVisitor() : vv(new std::vector<Vertex>()) {}

    void discover_vertex(Vertex v, const Graph &g)
    {
        if (boost::in_degree(v, g) != 0)
        {
            vv->push_back(v);
        }
        return;
    }
    std::vector<Vertex> &GetVector() const { return *vv; }

private:
    boost::shared_ptr<std::vector<Vertex>> vv;
};
// function declaration
bool Bellman_ford(Graph &g, std::vector<arc> edges, int number_of_nodes, weight_pmap weights, int *dist, E *pred, Vertex s);
inline void Update_pred(const Graph &g, Vertex v, bool *In_R, bool *reached_from_node_in_U, E *pred);

void check_sp(Graph &g, int number_of_nodes, Vertex s, weight_pmap weights, int *dist, E *pred);
void copy_graph(const LGraph &LG, Graph &BG, edge_array<int> &l_weight, WeightPrMap &b_weight,
                node_array<NodeInfo> &nodeInfo, DistanceMap &b_dist, VisitedMap &b_visit);
int NUM_OF_NODES = 5;
int NUM_OF_EDGES = 6;
int main()
{
    LGraph LG;
    edge_array<int> l_weight(LG);

    Graph BG;
    WeightPrMap b_weight = get(&EdgeProperties::weight, BG);
    DistanceMap b_dist = get(&NodeInfo::distance, BG);
    VisitedMap b_visit = get(&NodeInfo::visited, BG);

    random_graph(LG, NUM_OF_NODES, NUM_OF_EDGES, true, true, true);
    Make_Connected(LG);
    int i = 0;
    leda::edge e;

    l_weight.init(LG);
    leda::random_source S;
    std::vector<arc> edges_array;
    forall_edges(e, LG)
        l_weight[e] = (S(-100, -10));

    node_array<NodeInfo> nodeInfo(LG);
    node v;
    i = 0;
    forall_nodes(v, LG)
    {
        nodeInfo[v].visited = false;
        nodeInfo[v].distance = i++;
    }
    copy_graph(LG, BG, l_weight, b_weight, nodeInfo, b_dist, b_visit);

    std::pair<arcIt, arcIt> edge_it;
    int k = 0;
    Vertex src;
    Vertex trg;
    for (edge_it = edges(BG); edge_it.first != edge_it.second; ++edge_it.first)
    {
        src = source(*edge_it.first, BG);
        trg = target(*edge_it.first, BG);
        E temp_edge = E(src, trg);
        edges_array.push_back(*edge_it.first);
    }
    int dist[NUM_OF_NODES];
    for (i = 1; i < NUM_OF_NODES + 1; i++)
        dist[i] = INT_MAX;

    E pred[NUM_OF_NODES];
    Graph correctedGraph;

    // Vertex E = (Vertex)nullVertex;

    bool flag = true;
    flag = Bellman_ford(BG, edges_array, NUM_OF_NODES, b_weight, dist, pred, src);
    std::cout << "flag is: " << flag << std::endl;

    check_sp(BG, NUM_OF_NODES, src, b_weight, dist, pred);
    return 0;
}

bool Bellman_ford(Graph &g, std::vector<arc> edge_array, int number_of_nodes, weight_pmap weights, int *dist, E *pred, Vertex z)
{
    Vertex myNull = (Vertex)nullVertex;
    E myEdge = E(nullVertex, nullVertex);
    graph_traits<Graph>::vertex_iterator vi,
        vi_end;
    Vertex s = z;

    int phase_count = 0;
    std::queue<Vertex> Q;
    bool in_Q[number_of_nodes];
    int j = 0;

    // in_Q = false
    for (j = 0; j < number_of_nodes; j++)
        in_Q[j] = false;

    // pred = NULL
    for (boost::tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi, ++j)
        pred[*vi] = myEdge;

    Vertex temp;

    dist[s] = 0;
    Q.push(s);
    in_Q[s] = true;

    Q.push(myNull);

    while (phase_count < number_of_nodes)
    {

        temp = Q.front();
        Q.pop();
        if (temp == myNull)
        {
            phase_count++;
            if (Q.empty())
                return true;

            Q.push(myNull);
            continue;
        }
        else
        {
            in_Q[temp] = false;
        }
        int du = dist[temp];

        Vertex v, trg, src;

        std::pair<aIt, aIt> adj_edge;

        for (adj_edge = adjacent_vertices(temp, g); adj_edge.first != adj_edge.second; ++adj_edge.first)
        {
            std::pair<arc, bool> edge_temp = boost::edge(temp, *adj_edge.first, g);

            int d = du + weights[edge_temp.first];
            if ((pred[*adj_edge.first] == myEdge && *adj_edge.first != s) || (d < dist[*adj_edge.first])) //edw
            {
                dist[*adj_edge.first] = d;
                pred[*adj_edge.first] = E(temp, *adj_edge.first);
                if (!in_Q[*adj_edge.first])
                {
                    Q.push(*adj_edge.first);
                    in_Q[*adj_edge.first] = true;
                }
            }
        }
    }
    // BF: postprocessing

    if (pred[s] != myEdge)
        return false;

    bool in_R[number_of_nodes];
    for (j = 1; j < number_of_nodes + 1; j++)
        in_R[j] = false;

    graph_traits<Graph>::edge_iterator edge_it, edge_it_end;
    graph_traits<Graph>::vertex_descriptor u1, u2;
    graph_traits<Graph>::edge_descriptor e1;

    E temp_edge;
    std::vector<E> hidden_edge;
    int k = 0;

    for (int i = 0; i < edge_array.size(); i++)
    {
        Vertex src, trg;
        src = source(edge_array[i], g);
        trg = target(edge_array[i], g);
        temp_edge = E(src, trg);
        if (temp_edge != pred[trg])
        {
            hidden_edge.push_back(temp_edge);
            remove_edge(edge_array[i], g);
        }
    }

    //---------------CALL DFS----------------//
    MyVisitor vis;
    boost::depth_first_search(g, boost::visitor(vis));
    std::vector<Vertex> visited_vectrices = vis.GetVector();

    for (auto it = visited_vectrices.cbegin(); it != visited_vectrices.cend(); ++it)
        in_Q[*it] = true;
    // RESTORE EDGE
    for (int k = 0; k < hidden_edge.size(); k++)
    {
        add_edge(hidden_edge[k].first, hidden_edge[k].second, g);
    }

    bool reached_from_node_in_U[number_of_nodes];
    for (int i = 0; i < number_of_nodes; i++)
        reached_from_node_in_U[i] = false;

    k = 0;

    for (boost::tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi, ++k)
    {
        if (in_Q[*vi] && !reached_from_node_in_U[*vi])
            Update_pred(g, *vi, in_R, reached_from_node_in_U, pred);
    }

    return false;
}
inline void Update_pred(const Graph &g, Vertex v, bool *in_R, bool *reached_from_node_in_U, E *pred)
{
    reached_from_node_in_U[v] = true;

    Vertex trg;
    E temp_edge;

    std::pair<aIt, aIt> adj_edge;
    for (adj_edge = adjacent_vertices(v, g); adj_edge.first != adj_edge.second; ++adj_edge.first)
    {
        trg = *adj_edge.first;
        E temp_edge = E(v, trg);
        if (!reached_from_node_in_U[trg])
        {
            if (in_R[trg])
                pred[trg] = temp_edge;
            Update_pred(g, trg, in_R, reached_from_node_in_U, pred);
        }
    }
}

void check_sp(Graph &g, int number_of_nodes, Vertex s, weight_pmap weights, int *dist, E *pred)
{
    E myEdge = E(nullVertex, nullVertex);

    enum
    {
        NEG_CYCLE = -2,
        ATT_TO_CYCLE = -1,
        FINITE = 0,
        PLUS = 1,
        CYCLE = 2,
        ON_CUR_PATH = 3,
        UNKNOWN = 4
    };

    int label[number_of_nodes];
    bool reachable[number_of_nodes];

    for (int i = 0; i < number_of_nodes; i++)
    {
        reachable[i] = false;
        label[i] = UNKNOWN;
    }

    MyVisitor reached;
    boost::depth_first_search(g, boost::visitor(reached));

    std::vector<Vertex> visited_vectrices = reached.GetVector();
    graph_traits<Graph>::vertex_iterator vi, vi_end;

    for (std::vector<Vertex>::iterator it = visited_vectrices.begin(); it != visited_vectrices.end(); ++it)
        reachable[*it] = true;

    for (boost::tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
    {
        if (*vi != s)
        {
            assert((pred[*vi] == myEdge) == (reachable[*vi] == false));
            if (reachable[*vi] == false)
                label[*vi] = PLUS;
        }
    }
    // classification

    if (pred[s] == myEdge)
        label[s] = FINITE;
    for (boost::tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
    {
        if (label[*vi] == UNKNOWN)
        {
            std::queue<Vertex> S;
            Vertex w = *vi;
            while (label[w] == UNKNOWN)
            {
                label[w] = ON_CUR_PATH;
                S.push(w);
                w = source(pred[w], g);
            }
            // label all nodes on current path
            int t = label[w];
            if (t == ON_CUR_PATH)
            {
                Vertex z;
                do
                {
                    z = S.front();
                    S.pop();
                    label[z] = CYCLE;
                } while (z != w);
                while (!S.empty())
                {
                    label[S.front()] = ATT_TO_CYCLE;
                    S.pop();
                }
            }
            else
            {
                if (t == CYCLE)
                    t = ATT_TO_CYCLE;
                while (!S.empty())
                {
                    label[S.front()] = t;
                    S.pop();
                }
            }
        }
    }

    // condition 2
    for (boost::tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
    {
        if (label[*vi] == CYCLE)
        {
            Vertex w = *vi;
            int cycle_length = 0;
            do
            {
                Vertex trg = target(pred[w], g);
                Vertex src = source(pred[w], g);
                std::vector<E> hidden_edge;
                std::pair<arc, bool> e1 = boost::edge(src, trg, g);

                cycle_length += weights[e1.first];
                label[w] = NEG_CYCLE;
                // w = source(pred[w], g);
            } while (w != *vi);
            assert(cycle_length < 0);
        }
    }
    graph_traits<Graph>::edge_iterator edge_it, edge_it_end;
    graph_traits<Graph>::edge_descriptor e1;

    // condition 3, 4 and 5
    if (label[s] == FINITE)
        assert(dist[s] == 0);
    Vertex src, trg;
    for (boost::tie(edge_it, edge_it_end) = edges(g); edge_it != edge_it_end; ++edge_it)
    {
        src = source(*edge_it, g);
        trg = target(*edge_it, g);
        if (label[trg] == FINITE)
        {
            assert(label[src] == FINITE || label[src] == PLUS);
            if (label[src] == FINITE)
            {
                assert(dist[src] + weights[*edge_it] >= dist[trg]);
                if (E(src, trg) == pred[trg])
                    assert(dist[src] + weights[*edge_it] == dist[trg]);
            }
        }
    }
}
void copy_graph(const LGraph &LG, Graph &BG, edge_array<int> &l_weight, WeightPrMap &b_weight,
                node_array<NodeInfo> &nodeInfo, DistanceMap &b_dist, VisitedMap &b_visit)
{
    leda::node_array<Vertex> copy_in_BG(LG);
    arc a;
    leda::edge e;
    leda::node v;

    forall_nodes(v, LG)
    {
        copy_in_BG[v] = add_vertex(BG);
        b_dist[copy_in_BG[v]] = nodeInfo[v].distance;
        b_visit[copy_in_BG[v]] = nodeInfo[v].visited;
    }
    bool isAdded;
    forall_edges(e, LG)
    {
        tie(a, isAdded) = add_edge(copy_in_BG[source(e)], copy_in_BG[target(e)], BG);
        b_weight[a] = l_weight[e];
    }
}
