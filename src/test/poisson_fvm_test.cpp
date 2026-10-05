#include "cfd/fvm/fvm_dfdn.hpp"
#include "cfd/fvm/fvm_extended_collocations.hpp"
#include "cfd/fvm/fvm_gradient.hpp"
#include "cfd/grid/grid1d.hpp"
#include "cfd/grid/unstructured_grid2d.hpp"
#include "cfd/grid/vtk.hpp"
#include "cfd/mat/csrmat.hpp"
#include "cfd/mat/lodmat.hpp"
#include "cfd/mat/sparse_matrix_solver.hpp"
#include "cfd27_test.hpp"
#include "utils/filesystem.hpp"

using namespace cfd;

namespace {

class BaseProblem {
public:
    virtual ~BaseProblem() = default;
    virtual double exact_solution(Point p) const = 0;
    virtual double f(Point p) const = 0;

    // returns norm2(u - u_exact)
    double solve() {
        // 1. build SLAE
        CsrMatrix mat = approximate_lhs();
        std::vector<double> rhs = approximate_rhs();

        // 2. solve SLAE
        AmgcMatrixSolver solver; // iterative SLAE solver
        solver.set_matrix(mat);
        u_ = std::vector<double>(grid_->n_cells(), 0.0); // initial approximation for iterative solver
        solver.solve(rhs, u_);

        // 3. compute norm2
        return compute_norm2();
    }

    // saves numerical and exact solution into the vtk format
    void save_vtk(const std::string& filename) const {
        // save grid
        grid_->save_vtk(filename);

        // save numerical solution
        VtkUtils::add_cell_data(u_, "numerical", filename, grid_->n_cells());

        // save exact solution
        std::vector<double> exact(grid_->n_cells());
        for (size_t i = 0; i < grid_->n_cells(); ++i) {
            exact[i] = exact_solution(grid_->cell_center(i));
        }
        VtkUtils::add_cell_data(exact, "exact", filename);
    }

    std::shared_ptr<const IGrid> grid() const {
        return grid_;
    }

protected:
    BaseProblem(std::shared_ptr<IGrid> grid) : grid_(grid){};

    std::shared_ptr<IGrid> grid_;
    std::vector<double> u_;

    virtual CsrMatrix approximate_lhs() const {
        LodMatrix mat(grid_->n_cells());

        for (size_t icell = 0; icell < grid_->n_cells(); ++icell) {
            for (size_t iface: grid_->tab_cell_face(icell)) {
                auto [cell_left, cell_right] = grid_->tab_face_cell(iface);
                size_t jcell = (cell_left == icell) ? cell_right : cell_left;

                if (jcell != INVALID_INDEX) {
                    // icell and jcell exist => internal face
                    Point xi = grid_->cell_center(icell);
                    Point xj = grid_->cell_center(jcell);
                    double h = vector_abs(xj - xi);
                    double coef = grid_->face_area(iface) / h;
                    mat.add_value(icell, icell, coef);
                    mat.add_value(icell, jcell, -coef);
                } else {
                    // jcell doesn't exist => boundary face
                    Point gs = grid_->face_center(iface);
                    Point xi = grid_->cell_center(icell);
                    double h = vector_abs(gs - xi);
                    double coef = grid_->face_area(iface) / h;
                    mat.add_value(icell, icell, coef);
                }
            }
        }

        return mat.to_csr();
    }

    virtual std::vector<double> approximate_rhs() const {
        std::vector<double> rhs(grid_->n_cells(), 0.0);
        // internal
        for (size_t icell = 0; icell < grid_->n_cells(); ++icell) {
            double value = f(grid_->cell_center(icell));
            double volume = grid_->cell_volume(icell);
            rhs[icell] = value * volume;
        }
        // boundary faces
        for (size_t iface: grid_->boundary_faces()) {
            size_t icell = grid_->boundary_face_info(iface).icell;
            Point gs = grid_->face_center(iface);
            Point xi = grid_->cell_center(icell);
            double h = vector_abs(gs - xi);
            double coef = grid_->face_area(iface) / h;
            rhs[icell] += exact_solution(gs) * coef;
        }
        return rhs;
    }

    double compute_norm2() const {
        double norm2 = 0;
        double full_area = 0;
        for (size_t icell = 0; icell < grid_->n_cells(); ++icell) {
            double diff = u_[icell] - exact_solution(grid_->cell_center(icell));
            norm2 += grid_->cell_volume(icell) * diff * diff;
            full_area += grid_->cell_volume(icell);
        }
        return std::sqrt(norm2 / full_area);
    }
};

/////////////////////////////////////////////////
// 1D Worker
/////////////////////////////////////////////////
class Problem1D : public BaseProblem {
public:
    Problem1D(std::shared_ptr<IGrid1D> grid) : BaseProblem(grid) {}

    // u(x) = sin(10*x^2)
    double exact_solution(Point p) const override {
        double x = p.x;
        return sin(10 * x * x);
    }

    // -d^2 u(x)/d x^2
    double f(Point p) const override {
        double x = p.x;
        return 400 * x * x * sin(10 * x * x) - 20 * cos(10 * x * x);
    }
};

} // namespace

TEST_CASE("Poisson 1D solver, Finite Volume Method", "[poisson1-fvm]") {
    std::cout << std::endl << "--- [poisson1-fvm] --- " << std::endl;

    // precalculated norm2 results for some n_cells values
    // used for CHECK procedures
    std::map<size_t, double> norm2_for_compare{
        {10, 0.106539},
        {100, 0.00101714},
        {1000, 1.01641e-05},
    };

    // loop over n_cells value
    for (size_t n_cells: {10, 20, 50, 100, 200, 500, 1000}) {
        // build test solver
        auto grid = std::make_shared<Grid1D>(0, 1, n_cells);
        Problem1D worker(grid);

        // solve and find norm2
        double n2 = worker.solve();

        // save into poisson1_ncells={n_cells}.vtk
        worker.save_vtk("poisson1_fvm_n=" + std::to_string(n_cells) + ".vtk");

        // print (N_CELLS, NORM2) table entry
        std::cout << n_cells << " " << n2 << std::endl;

        // CHECK if result for this n_cells
        // presents in the norm2_for_compare dictionary
        auto found = norm2_for_compare.find(n_cells);
        if (found != norm2_for_compare.end()) {
            CHECK(n2 == Approx(found->second).margin(1e-6));
        }
    }
}
