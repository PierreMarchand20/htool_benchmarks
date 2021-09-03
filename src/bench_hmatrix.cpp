#include "bench_hmatrix.hpp"

int main(int argc, char *argv[]) {
    // Check the number of parameters
    if (argc < 1) {
        // Tell the user how to run the program
        std::cerr << "Usage: " << argv[0] << " n clustering type vectorisation compressor minclustersize outputpath" << std::endl;
        /* "Usage messages" are a conventional way of telling the user
        * how to run a program if they enter the command incorrectly.
        */
        return 1;
    }

    int n                    = std::stoi(argv[1]);
    std::string cluster_type = argv[2];
    char type                = *(argv[3]);         // S (float), D (double), C(complex<float>), Z (complex<double>)
    int vectorisation        = std::stoi(argv[4]); // 0 (no vectorisation), 1 (vector class library), 2 (xsimd unaligned), 2 (xsimd aligned)
    std::string compressor   = argv[5];
    int minclustersize       = std::stoi(argv[6]);
    std::string outputpath   = argv[7];

    // Initialize the MPI environment
    MPI_Init(&argc, &argv);
    int MPI_size;
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_size);

    // Create Geometries
    std::vector<double> p1(3 * n), p2(3 * n);
    double shift = 1e-5;
    double k     = create_geometry(n, p1.data());
    create_geometry(n, p2.data(), std::array<double, 3>{{shift, shift, shift}});
    if (type == 'S' || type == 'D') {
        k = 0;
    }

    // Clustering
    std::shared_ptr<htool::VirtualCluster> t, s;
    if (cluster_type == "PCARegularClustering") {
        t = std::make_shared<htool::Cluster<htool::PCARegularClustering>>();
        s = std::make_shared<htool::Cluster<htool::PCARegularClustering>>();
    } else if (cluster_type == "PCAGeometricClustering") {
        t = std::make_shared<htool::Cluster<htool::PCAGeometricClustering>>();
        s = std::make_shared<htool::Cluster<htool::PCAGeometricClustering>>();
    } else if (cluster_type == "BoundingBox1RegularClustering") {
        t = std::make_shared<htool::Cluster<htool::BoundingBox1RegularClustering>>();
        s = std::make_shared<htool::Cluster<htool::BoundingBox1RegularClustering>>();
    } else if (cluster_type == "BoundingBox1GeometricClustering") {
        t = std::make_shared<htool::Cluster<htool::BoundingBox1GeometricClustering>>();
        s = std::make_shared<htool::Cluster<htool::BoundingBox1GeometricClustering>>();
    } else {
        std::cerr << "Clustering not supported" << std::endl;
    }
    t->set_minclustersize(minclustersize);
    s->set_minclustersize(minclustersize);
    t->build(n, p1.data());
    s->build(n, p2.data());

    std::tuple<double, double, double, double, double> results;
    switch (type) {
    case 'S':
        if (compressor == "SVD") {
            results = bench_hmatrix<float, htool::SVD>(n, k, p1, p2, t, s, vectorisation);
        } else if (compressor == "fullACA") {
            results = bench_hmatrix<float, htool::fullACA>(n, k, p1, p2, t, s, vectorisation);
        } else if (compressor == "partialACA") {
            results = bench_hmatrix<float, htool::partialACA>(n, k, p1, p2, t, s, vectorisation);
        } else if (compressor == "sympartialACA") {
            results = bench_hmatrix<float, htool::sympartialACA>(n, k, p1, p2, t, s, vectorisation);
        }
        break;
    case 'D':
        if (compressor == "SVD") {
            results = bench_hmatrix<double, htool::SVD>(n, k, p1, p2, t, s, vectorisation);
        } else if (compressor == "fullACA") {
            results = bench_hmatrix<double, htool::fullACA>(n, k, p1, p2, t, s, vectorisation);
        } else if (compressor == "partialACA") {
            results = bench_hmatrix<double, htool::partialACA>(n, k, p1, p2, t, s, vectorisation);
        } else if (compressor == "sympartialACA") {
            results = bench_hmatrix<double, htool::sympartialACA>(n, k, p1, p2, t, s, vectorisation);
        }

        break;
    case 'C':
        if (compressor == "SVD") {
            results = bench_hmatrix<std::complex<float>, htool::SVD>(n, k, p1, p2, t, s, vectorisation);
        } else if (compressor == "fullACA") {
            results = bench_hmatrix<std::complex<float>, htool::fullACA>(n, k, p1, p2, t, s, vectorisation);
        } else if (compressor == "partialACA") {
            results = bench_hmatrix<std::complex<float>, htool::partialACA>(n, k, p1, p2, t, s, vectorisation);
        } else if (compressor == "sympartialACA") {
            results = bench_hmatrix<std::complex<float>, htool::sympartialACA>(n, k, p1, p2, t, s, vectorisation);
        }

        break;
    case 'Z':
        if (compressor == "SVD") {
            results = bench_hmatrix<std::complex<double>, htool::SVD>(n, k, p1, p2, t, s, vectorisation);
        } else if (compressor == "fullACA") {
            results = bench_hmatrix<std::complex<double>, htool::fullACA>(n, k, p1, p2, t, s, vectorisation);
        } else if (compressor == "partialACA") {
            results = bench_hmatrix<std::complex<double>, htool::partialACA>(n, k, p1, p2, t, s, vectorisation);
        } else if (compressor == "sympartialACA") {
            results = bench_hmatrix<std::complex<double>, htool::sympartialACA>(n, k, p1, p2, t, s, vectorisation);
        }

        break;
    default:
        std::cerr << "type not supported" << std::endl;
        break;
    }

    // Save
    std::ifstream infile((outputpath + "bench_hmatrix.csv").c_str());
    std::ofstream output((outputpath + "bench_hmatrix.csv").c_str(), std::ios_base::app);
    if (!infile.good()) {
        output << "type,freq,mpi,threads,n,time assemble,time prod,space saving,compressor,cluster_type,vectorized,checksum" << std::endl;
    }
    output << type << "," << k << "," << MPI_size << "," << omp_thread_count() << "," << n << "," << std::get<0>(results) << "," << std::get<1>(results) << "," << std::get<4>(results) << "," << compressor << "," << cluster_type << "," << vectorisation_name(vectorisation) << ",";

    if (type == 'S' || type == 'D') {
        output << std::get<2>(results) << std::endl;
    } else if (type == 'C' || type == 'Z') {
        output << std::get<2>(results) << "+i" << std::get<3>(results) << std::endl;
    }

    // Finalize the MPI environment.
    MPI_Finalize();
}