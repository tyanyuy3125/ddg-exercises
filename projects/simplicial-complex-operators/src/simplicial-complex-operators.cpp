// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

#include <Eigen/Sparse>
#include <algorithm>
#include <set>

using namespace geometrycentral;
using namespace geometrycentral::surface;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v];
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {

    // TODO
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.
    Eigen::SparseMatrix<size_t> sparse_mat(mesh->nEdges(), mesh->nVertices());
    for (Edge e : mesh->edges()) {
        sparse_mat.insert(e.getIndex(), e.firstVertex().getIndex()) = 1;
        sparse_mat.insert(e.getIndex(), e.secondVertex().getIndex()) = 1;
    }

    return sparse_mat;
    // return identityMatrix<size_t>(1); // placeholder
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {

    // TODO
    Eigen::SparseMatrix<size_t> sparse_mat(mesh->nFaces(), mesh->nEdges());
    // return identityMatrix<size_t>(1); // placeholder
    // for(Face f: mesh->faces())
    // {
    //     for(Vertex v : f.adjacentVertices())
    //     {
    //         sparse_mat.insert(f.getIndex(), v.getIndex()) = 1;
    //     }
    // }
    for (Edge e : mesh->edges()) {
        for (Face f : e.adjacentFaces()) {
            sparse_mat.insert(f.getIndex(), e.getIndex()) = 1;
        }
    }
    return sparse_mat;
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.isComplex
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {

    // TODO
    // return Vector<size_t>::Zero(1);
    Vector<size_t> ret = Vector<size_t>::Zero(mesh->nVertices());
    for (auto vert : subset.vertices) {
    }
    return ret;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {

    // TODO
    // return Vector<size_t>::Zero(1);
    Vector<size_t> ret = Vector<size_t>::Zero(mesh->nEdges());
    return ret;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {

    // TODO
    // return Vector<size_t>::Zero(1);
    Vector<size_t> ret = Vector<size_t>::Zero(mesh->nFaces());
    return ret;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {

    // TODO
    auto vertices = subset.vertices;
    auto edges = subset.edges;
    auto expanded_vertices = vertices;
    auto expanded_edges = edges;

    // To be expanded and returned.
    auto faces = subset.faces;

    for (auto vertex : vertices) {
        auto vert_obj = mesh->vertex(vertex);
        for (auto adj_edge : vert_obj.adjacentEdges()) {
            expanded_edges.insert(adj_edge.getIndex());
        }
    }

    for (auto edge : expanded_edges) {
        auto edge_obj = mesh->edge(edge);
        for (auto adj_face : edge_obj.adjacentFaces()) {
            faces.insert(adj_face.getIndex());
        }
    }

    MeshSubset ret(expanded_vertices, expanded_edges, faces);
    return ret;
    // return subset; // placeholder
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {

    // TODO
    auto vertices = subset.vertices;
    auto edges = subset.edges;
    auto expanded_vertices = vertices;
    auto expanded_edges = edges;

    // To be expanded and returned.
    auto faces = subset.faces;

    for (auto face : faces) {
        auto face_obj = mesh->face(face);
        for (auto clo_edge : face_obj.adjacentEdges()) {
            expanded_edges.insert(clo_edge.getIndex());
        }
    }

    for (auto edge : expanded_edges) {
        auto edge_obj = mesh->edge(edge);
        for (auto clo_vert : edge_obj.adjacentVertices()) {
            expanded_vertices.insert(clo_vert.getIndex());
        }
    }

    MeshSubset ret(expanded_vertices, expanded_edges, faces);
    return ret;
    // return subset; // placeholder
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {

    // TODO
    // return subset; // placeholder
    MeshSubset star_set = closure(star(subset));
    MeshSubset closure_set = star(closure(subset));
    MeshSubset ret;
    auto set_difference = [](const std::set<size_t>& a, const std::set<size_t>& b, std::set<size_t>& target) -> void {
        std::set_difference(a.begin(), a.end(), b.begin(), b.end(), std::inserter(target, target.begin()));
    };
    set_difference(star_set.vertices, closure_set.vertices, ret.vertices);
    set_difference(star_set.edges, closure_set.edges, ret.edges);
    set_difference(star_set.faces, closure_set.faces, ret.faces);
    // return closure(ret);
    return ret;
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    // TODO
    MeshSubset clo_test = closure(subset);
    return (clo_test.vertices == subset.vertices) && (clo_test.edges == subset.edges) &&
           (clo_test.faces == subset.faces); // placeholder
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {

    // TODO
    if (!isComplex(subset)) {
        return -1;
    }
    if (subset.edges.empty() && subset.faces.empty()) {
        return 0;
    }
    if (subset.faces.empty()) {
        for (auto test_vertex : subset.vertices) {
            auto test_vertex_obj = mesh->vertex(test_vertex);
            // test_vertex_obj.adjacentEdges()
            for (auto test_vertex_adj_edge : test_vertex_obj.adjacentEdges()) {
                if(subset.edges.find(test_vertex_adj_edge.getIndex()) == subset.edges.end())
                {
                    return -1;
                }
            }
        }
        return 1;
    }
    return -1; // placeholder
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {

    // TODO
    return subset; // placeholder
}