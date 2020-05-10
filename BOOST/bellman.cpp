#include <boost/config.hpp>
#include <vector>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/lockfree/queue.hpp>
#include <queue>
#include <utility>
#include <boost/graph/bellman_ford_shortest_paths.hpp>

using namespace boost;

struct EdgeProperties
{
    int weight;
};

typedef adjacency_list<vecS, vecS, bidirectionalS, no_property, EdgeProperties> Graph;
typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef property_map<Graph, int EdgeProperties::*>::type weight_pmap;
typedef std::pair<int, int> E;
bool Bellman_ford(const Graph &g, int number_of_nodes, weight_pmap weights, int *dist, E *pred, Vertex s);
inline void Update_pred(const Graph &g, Vertex v, bool *In_R, bool *reached_from_node_in_U, E *pred);
int main()
{
    enum
    {
        u,
        v,
        x,
        y,
        z,
        N
    };
    char name[] = {'u', 'v', 'x', 'y', 'z'};
    const int n_edges = 10;
    E pred[N];
    E edge_array[] = {E(u, y), E(u, v), E(u, x), E(v, u), E(x, y), E(x, v), E(y, v), E(y, z), E(z, u), E(z, x)};
    int weight[n_edges] = {-4, 8, 5, -2, 9, -3, 7, 2, 6, 7};

    Graph g(edge_array, edge_array + n_edges, N); // Graph
    weight_pmap weights = get(&EdgeProperties::weight, g);
    graph_traits<Graph>::edge_iterator ei, ei_end;
    int i = 0;
    for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei, ++i)
        weights[*ei] = weight[i];

    typedef property_map<Graph, vertex_index_t>::type IndexMap;
    IndexMap index = get(vertex_index, g);
    std::vector<std::size_t> parent(N);
    int distance[N];

    for (i = 0; i < N; i++)
    {
        parent[i] = i;
        distance[i] = INT_MAX;
        pred[i] = {};
    }

    distance[z] = 0;
    Vertex s;
    Bellman_ford(g, N, weights, distance, pred, z);
    return 0;
}

bool Bellman_ford(const Graph &g, int number_of_nodes, weight_pmap weights, int *dist, E *pred, Vertex z)
{
    graph_traits<Graph>::vertex_iterator vi, vi_end;
    Vertex s = z;
    int phase_count = 0;

    std::queue<Vertex> Q;
    bool in_Q[number_of_nodes];
    int j = 0;
    for (j = 0; j < number_of_nodes; j++)
        in_Q[j] = false;
    for (boost::tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi, ++j)
        pred[*vi] = E(NULL, NULL);

    dist[s] = 0;
    Q.push(s);
    in_Q[s] = true;
    Q.push((Vertex)NULL);
    Vertex u;
    while (phase_count < number_of_nodes)
    {
        u = Q.front();
        Q.pop();
        if (u == 0)
        {
            phase_count++;
            if (Q.empty())
                continue;
            Q.push((Vertex)NULL);
            continue;
        }
        else
            in_Q[u] = false;
        int du = dist[u];
        boost::graph_traits<Graph>::in_edge_iterator in_edge, in_edge_end;
        Vertex v, src;
        for (boost::tie(in_edge, in_edge_end) = in_edges(u, g); in_edge != in_edge_end; ++in_edge)
        {
            v = target(*in_edge, g);
            src = source(*in_edge, g);
            int d = du + weights[*in_edge];
            if ((pred[v] == E(NULL, NULL) && v != s) || (d < dist[v])) //edw
            {
                dist[v] = d;
                pred[v] = E(src, v);
                if (!in_Q[v])
                {
                    Q.push(v);
                    in_Q[v] = true;
                }
            }
        }
    }
    if (pred[s] != E(NULL, NULL))
        return false;
    bool in_R[number_of_nodes];
    for (j = 0; j < number_of_nodes; j++)
        in_R[j] = false;
    graph_traits<Graph>::edge_iterator edge_it, edge_it_end;
    int k = 0;
    Vertex u1, u2;
    E temp_edge;
    for (boost::tie(edge_it, edge_it_end) = edges(g); edge_it != edge_it_end; ++edge_it, ++k)
    {
        u1 = source(*edge_it, g);
        u2 = target(*edge_it, g);
        temp_edge = E(u1, u2);
        if (temp_edge != pred[target(*edge_it, g)])
        {
            // hide edge???
        }
        //CALL DFS

        // RESTORE EDGE

        bool reached_from_node_in_U[number_of_nodes];
        graph_traits<Graph>::vertex_iterator vi, vi_end;
        k = 0;
        for (boost::tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi, ++k)
        {
            if (in_Q[*vi] && !reached_from_node_in_U[*vi])
                Update_pred(g, *vi, in_R, reached_from_node_in_U, pred);
        }
    }
    return false;
}
inline void Update_pred(const Graph &g, Vertex v, bool *in_R, bool *reached_from_node_in_U, E *pred)
{
    reached_from_node_in_U[v] = true;
    boost::graph_traits<Graph>::in_edge_iterator in_edge, in_edge_end;

    Vertex src, trg;
    E temp_edge;

    for (boost::tie(in_edge, in_edge_end) = in_edges(v, g); in_edge != in_edge_end; ++in_edge)
    {
        trg = target(*in_edge, g);
        src = source(*in_edge, g);
        temp_edge = E(src, trg);

        if (!reached_from_node_in_U[trg])
        {
            if (in_R[trg])
                pred[trg] = temp_edge;
            Update_pred(g, trg, in_R, reached_from_node_in_U, pred);
        }
    }
}
