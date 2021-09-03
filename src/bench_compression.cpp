#include "bench_compression.hpp"

int main(int argc, char *argv[]) {

    // Check the number of parameters
    if (argc < 1) {
        // Tell the user how to run the program
        std::cerr << "Usage: " << argv[0] << " n m clustering type vectorisation compressor outputpath" << std::endl;
        /* "Usage messages" are a conventional way of telling the user
        * how to run a program if they enter the command incorrectly.
        */
        return 1;
    }

    int nr                   = std::stoi(argv[1]);
    int nc                   = std::stoi(argv[2]);
    std::string cluster_type = argv[3];
    char type                = *(argv[4]);
    int vectorisation        = std::stoi(argv[5]);
    std::string compressor   = argv[6];
    std::string outputpath   = argv[7];

    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    // Create Geometries
    std::vector<double> p1(3 * nr), p2(3 * nc);
    double shift = 10;
    double k     = create_geometry(nr, p1.data());
    create_geometry(nc, p2.data(), std::array<double, 3>{{shift, shift, shift}});
    if (type == 'S' || type == 'D') {
        k = 0;
    }

    // Clustering
    std::unique_ptr<htool::VirtualCluster> t, s;
    if (cluster_type == "PCARegularClustering") {
        t = std::make_unique<htool::Cluster<htool::PCARegularClustering>>();
        s = std::make_unique<htool::Cluster<htool::PCARegularClustering>>();
    } else if (cluster_type == "PCAGeometricClustering") {
        t = std::make_unique<htool::Cluster<htool::PCAGeometricClustering>>();
        s = std::make_unique<htool::Cluster<htool::PCAGeometricClustering>>();
    } else if (cluster_type == "BoundingBox1RegularClustering") {
        t = std::make_unique<htool::Cluster<htool::BoundingBox1RegularClustering>>();
        s = std::make_unique<htool::Cluster<htool::BoundingBox1RegularClustering>>();
    } else if (cluster_type == "BoundingBox1GeometricClustering") {
        t = std::make_unique<htool::Cluster<htool::BoundingBox1GeometricClustering>>();
        s = std::make_unique<htool::Cluster<htool::BoundingBox1GeometricClustering>>();
    } else {
        std::cerr << "Clustering not supported" << std::endl;
    }

    t->build(nr, p1.data());
    s->build(nc, p2.data());
    htool::VirtualCluster *cluster_target = t.get();
    htool::VirtualCluster *cluster_source = s.get();

    // Bench
    std::tuple<double, double, double, double, int> results;
    switch (type) {
    case 'S':
        if (compressor == "SVD") {
            results = bench_compression<float, htool::SVD>(nr, nc, k, p1, p2, cluster_target, cluster_source, vectorisation);
        } else if (compressor == "fullACA") {
            results = bench_compression<float, htool::fullACA>(nr, nc, k, p1, p2, cluster_target, cluster_source, vectorisation);
        } else if (compressor == "partialACA") {
            results = bench_compression<float, htool::partialACA>(nr, nc, k, p1, p2, cluster_target, cluster_source, vectorisation);
        } else if (compressor == "sympartialACA") {
            results = bench_compression<float, htool::sympartialACA>(nr, nc, k, p1, p2, cluster_target, cluster_source, vectorisation);
        }
        break;
    case 'D':
        if (compressor == "SVD") {
            results = bench_compression<double, htool::SVD>(nr, nc, k, p1, p2, cluster_target, cluster_source, vectorisation);
        } else if (compressor == "fullACA") {
            results = bench_compression<double, htool::fullACA>(nr, nc, k, p1, p2, cluster_target, cluster_source, vectorisation);
        } else if (compressor == "partialACA") {
            results = bench_compression<double, htool::partialACA>(nr, nc, k, p1, p2, cluster_target, cluster_source, vectorisation);
        } else if (compressor == "sympartialACA") {
            results = bench_compression<double, htool::sympartialACA>(nr, nc, k, p1, p2, cluster_target, cluster_source, vectorisation);
        }

        break;
    case 'C':
        if (compressor == "SVD") {
            results = bench_compression<std::complex<float>, htool::SVD>(nr, nc, k, p1, p2, cluster_target, cluster_source, vectorisation);
        } else if (compressor == "fullACA") {
            results = bench_compression<std::complex<float>, htool::fullACA>(nr, nc, k, p1, p2, cluster_target, cluster_source, vectorisation);
        } else if (compressor == "partialACA") {
            results = bench_compression<std::complex<float>, htool::partialACA>(nr, nc, k, p1, p2, cluster_target, cluster_source, vectorisation);
        } else if (compressor == "sympartialACA") {
            results = bench_compression<std::complex<float>, htool::sympartialACA>(nr, nc, k, p1, p2, cluster_target, cluster_source, vectorisation);
        }

        break;
    case 'Z':
        if (compressor == "SVD") {
            results = bench_compression<std::complex<double>, htool::SVD>(nr, nc, k, p1, p2, cluster_target, cluster_source, vectorisation);
        } else if (compressor == "fullACA") {
            results = bench_compression<std::complex<double>, htool::fullACA>(nr, nc, k, p1, p2, cluster_target, cluster_source, vectorisation);
        } else if (compressor == "partialACA") {
            results = bench_compression<std::complex<double>, htool::partialACA>(nr, nc, k, p1, p2, cluster_target, cluster_source, vectorisation);
        } else if (compressor == "sympartialACA") {
            results = bench_compression<std::complex<double>, htool::sympartialACA>(nr, nc, k, p1, p2, cluster_target, cluster_source, vectorisation);
        }

        break;
    default:
        std::cerr << "type not supported" << std::endl;
        break;
    }

    // Save
    std::ifstream infile((outputpath + "bench_compression.csv").c_str());
    std::ofstream output((outputpath + "bench_compression.csv").c_str(), std::ios_base::app);
    if (!infile.good()) {
        output << "type,freq,m,n,time assemble,time prod,rank,compressor,clustering,vectorized,checksum" << std::endl;
    }
    output << type << "," << k << "," << nr << "," << nc << "," << std::get<0>(results) << "," << std::get<1>(results) << "," << std::get<4>(results) << "," << compressor << "," << cluster_type << "," << vectorisation_name(vectorisation) << ",";

    if (type == 'S' || type == 'D') {
        output << std::get<2>(results) << std::endl;
    } else if (type == 'C' || type == 'Z') {
        output << std::get<2>(results) << "+i" << std::get<3>(results) << std::endl;
    }

    // Finalize the MPI environment.
    MPI_Finalize();
}