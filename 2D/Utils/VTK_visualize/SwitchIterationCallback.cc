/* ********************************************
 * SwitchIterationCallback class implementation
 * ******************************************** */
#include <string>
#include <sstream>

#include <vtkObject.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyle.h>

#include "include/Declare_functions.hh"
#include "include/SwitchIterationCallback.hh"

using namespace std;


/* =========================================================
 * Callback used to switch iteration on left/right key press
 * ========================================================= */
void SwitchIterationCallback::Execute(vtkObject *caller, unsigned long, void*) {
    const auto niters_m1 = this->niters - 1;
          auto istyle    = static_cast<vtkInteractorStyle*>(caller);
    const auto key       = static_cast<string>(istyle->GetInteractor()->GetKeySym());

    if (key == "Right") {
        if (this->current_iter == niters_m1) {
            this->current_iter = 0;
        } else if (this->current_iter < niters_m1) {
            ++this->current_iter;
        }
    } else if (key == "Left") {
        if (this->current_iter == 0) {
            this->current_iter = niters_m1;
        } else if (this->current_iter > 0) {
            --this->current_iter;
        }
    } else {
        return;
    }

    ostringstream text_ss;
    text_ss << "Iteration " << current_iter       << endl
            << "Left/Right: switch iteration"     << endl
            << "<Q>/<E>: quit"                    << endl
            << "Left mouse button + drag: scroll" << endl
            << "Right mouse button + drag or mouse wheel rolling: zoom in/out";
    this->text_actor->SetInput(text_ss.str().c_str());

    load_iteration(current_iter, this->nx, this->ny, this->file_id,
                   this->grid, this->warp, this->window);

    return;
}
