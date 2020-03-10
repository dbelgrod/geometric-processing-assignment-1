#include <iostream>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
/*** insert any libigl headers here ***/
#include <igl/vertex_triangle_adjacency.h>
#include <igl/adjacency_list.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_corner_normals.h>
#include <igl/facet_components.h>
#include <igl/jet.h>
#include <igl/edge_topology.h>
#include <igl/barycenter.h>
#include <math.h>
#include <igl/matlab_format.h>

#define PI 3.14159265

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;

// Vertex array, #V x3
Eigen::MatrixXd V;
// Face array, #F x3
Eigen::MatrixXi F;
// Per-face normal array, #F x3
Eigen::MatrixXd FN;
// Per-vertex normal array, #V x3
Eigen::MatrixXd VN;
// Per-corner normal array, (3#F) x3
Eigen::MatrixXd CN;
// Vectors of indices for adjacency relations
std::vector<std::vector<int> > VF, VFi, VV;
// Integer vector of component IDs per face, #F x1
Eigen::VectorXi cid;
// Per-face color array, #F x3
Eigen::MatrixXd component_colors_per_face;

template <typename IndexVector>
void print_vector(std::vector<std::vector<IndexVector> >& A){
    for (int i=0; i< A.size(); i++){
            for (int j=0; j<A[i].size(); j++){
                cout << A[i][j] << " ";
            }
            cout << endl;
            }
}


bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers) {
    if (key == '1') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing vertex to face relations here;
        // store in VF,VFi.
        igl::vertex_triangle_adjacency(V,F, VF, VFi);
        print_vector(VF);
        cout << "Printing vertex to face for #F:" << F.rows() << endl;
    }

    if (key == '2') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing vertex to vertex relations here:
        // store in VV.
        igl::adjacency_list(F, VV);
        print_vector(VV);
        cout << "Printing vertex adjancency for #V:" << V.rows() 
                          << "==" << VV.size() << endl;
    }

    if (key == '3') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        FN.setZero(F.rows(),3);
        // Add your code for computing per-face normals here: store in FN.
        igl::per_face_normals(V,F, FN);
        // Set the viewer normals.
        viewer.data().set_normals(FN);
    }

    if (key == '4') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing per-vertex normals here: store in VN.
        igl::per_vertex_normals(V,F, VN);
        // Set the viewer normals.
        viewer.data().set_normals(VN);
    }

    if (key == '5') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing per-corner normals here: store in CN.
        igl::per_corner_normals(V,F, 120, CN);
        //Set the viewer normals
        viewer.data().set_normals(CN);
    }

    if (key == '6') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        component_colors_per_face.setZero(F.rows(),3);
        // Add your code for computing per-face connected components here:
        // store the component labels in cid.
        
        igl::facet_components(F, cid);
        //igl::jet

        // Compute colors for the faces based on components, storing them in
        // component_colors_per_face.
        igl::jet(cid, true, component_colors_per_face);
        // Set the viewer colors
        viewer.data().set_colors(component_colors_per_face);
        cout << cid << endl;
    }

    if (key == '7') {
        Eigen::MatrixXd Vout=V;
        Eigen::MatrixXi Fout=F;
        // Add your code for sqrt(3) subdivision here.
        
        Eigen::MatrixXd M; //(F.rows(), V.cols());
        igl::barycenter(V,F,M);

        vector<vector<int> > adj;
        igl::adjacency_list(F, adj);

        Eigen::MatrixXd P;
        P = Eigen::MatrixXd::Constant(V.rows(), V.cols(), 0);

        for (int i=0; i<adj.size(); i++){ //for reach  vertex
            int n = adj[i].size(); //n is size of each list
            float a = (4 - cos(2*PI / n)) / 9;

            for (int j=0; j< adj[i].size(); j++){ //for each list
                for (int k=0; k< P.cols(); k++){ //for every col in P
                P(i, k) += V(adj[i][j],k);
                }
            }
            for (int k=0; k<P.cols();k++){
                P(i,k) = P(i,k) * a / n + (1 - a) * V(i,k);
            }
        }

        Eigen::MatrixXi EV;
        Eigen::MatrixXi FE;
        Eigen::MatrixXi EF;
        igl::edge_topology(P, F, EV, FE, EF); //switched from V to P
        Eigen::MatrixXi Ffin(3*F.rows(), F.cols());
        Eigen::MatrixXd Vfin(P.rows() + M.rows(), P.cols());
        Vfin << P,   
                M;
        
        //std::cout<<igl::matlab_format(Vfin,"Vfin")<< std::endl;
        //std::cout<<igl::matlab_format(V,"V")<< std::endl;

        //std::vector<int> seen_faces;
        int k = 0;
        for (int i=0; i<EF.rows(); i++){
            //we get EV which the 2 vertices that make up the edge
            int v1 = EV(i,0);
            int v2 = EV(i,1); 

            //must be 2 faces only per edge?
            int f1 = EF(i,0);
            int f2 = EF(i,1);

            //we can just take the face index then go find in Vfin

            Ffin(k,0) = v1;
            Ffin(k,1) = V.rows()+f2;
            Ffin(k,2) = V.rows()+f1;
            k++;

            Ffin(k,0) = v2;
            
            Ffin(k,1) = V.rows()+f1;
            Ffin(k,2) = V.rows()+f2; //fix orientation
            k++;
            
        }
       
        Fout = Ffin; Vout = Vfin;
        std::cout<<igl::matlab_format(Fout,"F")<< std::endl;
        // Set up the viewer to display the new mesh
        V = Vout; F = Fout;
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
    }

    return true;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        cout << "Usage ex1_bin mesh.obj" << endl;
        exit(0);
    }

    // Read mesh
    igl::readOFF(argv[1],V,F);

    // Show the mesh
    Viewer viewer;
    viewer.callback_key_down = callback_key_down;
    
    
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);
    
    menu.callback_draw_viewer_menu = [&](){
        menu.draw_viewer_menu();
    };

    callback_key_down(viewer, '1', 0);

    viewer.launch();
}
