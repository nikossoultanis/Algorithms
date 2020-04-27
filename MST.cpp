#include <LEDA/graph/graph.h>
#include <LEDA/core/random_source.h>
#include <LEDA/core/p_queue.h>
#include <LEDA/core/array.h>
#include "LEDA/system/timer.h"
#include <LEDA/core/dynamic_trees.h>
#include <LEDA/graph/min_span.h>

#include <list>
#include <math.h>
#define NUM_NODES 16000

using namespace leda;

list<edge> MST_finder(const graph &G, const edge_array<int> &cost);
bool MST_checker(const graph &G, const list<edge> &T, const edge_array<int> &cost, const list<edge> &allEdges);

struct myNode
{
    int size;
    node first;
    node last;
    std::list<node> myNodeList;
};

int main()
{
    graph G;
    random_source S;
    random_graph(G, NUM_NODES, 2 * NUM_NODES * log10(NUM_NODES), true, true, true);
    Make_Connected(G);
    // complete_ugraph(G, 3500);
    G.make_undirected();
    list<edge> allEdges; // total edges
    edge e;
    edge_array<int> cost(G);
    forall_edges(e, G)
    {
        cost[e] = (S(10, 10000));
        allEdges.append(e);
    }

    bool flag;
    list<edge> T;
    timer dfs;
    dfs.start();
    T = MST_finder(G, cost);
    dfs.stop();
    std::cout << "MST_finder() time " << dfs.elapsed_time() << "\n";

    dfs.start();
    MIN_SPANNING_TREE(G, cost);
    dfs.stop();
    std::cout << "LEDA MIN_SPANNING_TREE() time " << dfs.elapsed_time() << "\n";

    // flag = MST_checker(G, T, cost, allEdges);
    // std::cout << flag << std::endl;
    return 0;
}

list<edge> MST_finder(const graph &G, const edge_array<int> &cost)
{
    list<edge> T;

    node_array<myNode> allNodes(G);
    p_queue<int, edge> PQ;
    edge e;
    node n;
    forall_edges(e, G)
    {
        PQ.insert(cost[e], e);
    }
    forall_nodes(n, G)
    {
        myNode obj;
        obj.size = 1;
        obj.first = n;
        obj.last = n;
        obj.myNodeList.push_back(n);
        allNodes[n] = obj;
    }
    while (!PQ.empty())
    {
        e = PQ.inf(PQ.find_min());
        PQ.del_min();

        node src = G.source(e);
        node trg = G.target(e);
        if (allNodes[src].first != allNodes[trg].first)
        {
            T.append(e);
            if (allNodes[src].size >= allNodes[trg].size)
            {
                node temp = allNodes[trg].first;

                forall_nodes(n, G)
                {
                    if (allNodes[n].first == temp)
                    {
                        allNodes[n].first = allNodes[src].first;
                        // allNodes[src].myNodeList.push_back(n);
                        allNodes[src].size = allNodes[src].size + allNodes[trg].size;
                        allNodes[trg].size = allNodes[src].size;
                        // allNodes[src].last = allNodes[trg].last;
                    }
                }
            }
            else
            {
                node temp = allNodes[src].first;
                forall_nodes(n, G)
                {
                    if (allNodes[n].first == temp)
                    {
                        allNodes[n].first = allNodes[trg].first;
                        // allNodes[src].myNodeList.push_back(n);
                        allNodes[trg].size = allNodes[trg].size + allNodes[src].size;
                        allNodes[src].size = allNodes[trg].size;
                        // allNodes[trg].last = allNodes[src].last;
                    }
                }
            }
        }
    }
    std::cout << "MST has: " << T.size() << " edges, as the number of nodes - 1" << std::endl;

    return T;
}

bool MST_checker(const graph &G, const list<edge> &T, const edge_array<int> &cost, const list<edge> &allEdges)
{
    list<edge> nonTreeEdges; // not tree edges
    node_array<vertex> allTrees(G);
    dynamic_trees D;

    int counter = 0;

    edge tree_edge;
    edge e;
    // find all non tree edges
    forall(e, allEdges)
    {
        forall(tree_edge, T)
        {
            if (e == tree_edge)
                counter++;
        }

        if (counter != 1)
        {
            nonTreeEdges.append(e);
        }
        counter = 0;
    }
    node n;
    // vertex for each node
    forall_nodes(n, G)
    {
        allTrees[n] = D.make();
    }
    node src;
    node trg;
    // linking
    forall(e, T)
    {
        src = G.source(e);
        trg = G.target(e);
        if (D.root(allTrees[src]) == allTrees[src])
            D.link(allTrees[src], allTrees[trg], cost[e]);
        else
        {
            D.evert(allTrees[src]);
            D.link(allTrees[src], allTrees[trg], cost[e]);
        }
    }

    bool flag = true;
    // searching through the tree
    forall(e, nonTreeEdges)
    {
        src = G.source(e);
        trg = G.target(e);
        vertex temp_src = allTrees[src];
        vertex temp_trg = allTrees[trg];
        vertex common_ancestor;
        common_ancestor = D.lca(temp_src, temp_trg);
        while (temp_src != common_ancestor && flag != false)
        {

            if (D.cost(temp_src) <= cost[e]) // from source climb higher
            {
                temp_src = D.parent(temp_src); // climb to the parent. Stop if you find common ancestor
            }
            else
            {
                flag = false;
                return flag;
            }
        }
        while (temp_trg != common_ancestor != 0 && flag != false) // from target climb higher
        {
            if (D.cost(temp_trg) <= cost[e])
            {
                temp_trg = D.parent(temp_trg); // climb to the parent.Stop if you find common ancestor
            }
            else
            {
                flag = false;
                return flag;
            }
        }
    }
    return flag;
}