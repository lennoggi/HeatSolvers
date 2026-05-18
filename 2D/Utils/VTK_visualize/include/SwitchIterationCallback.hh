#ifndef SWITCH_ITERATION_CALLBACK
#define SWITCH_ITERATION_CALLBACK

#include <hdf5.h>

#include <vtkObject.h>
#include <vtkCommand.h>
#include <vtkRectilinearGrid.h>
#include <vtkWarpScalar.h>
#include <vtkRenderWindow.h>


class SwitchIterationCallback : public vtkCommand {
    public:
        hid_t               file_id      = -1;
        size_t              current_iter = 0;
        size_t              niters       = 0;
        size_t              nx           = 0;
        size_t              ny           = 0;
        vtkRectilinearGrid* grid         = nullptr;
        vtkWarpScalar*      warp         = nullptr;
        vtkRenderWindow*    window       = nullptr;

        static SwitchIterationCallback* New() {
            return new SwitchIterationCallback;
        }

        void Execute(vtkObject *caller, unsigned long, void*) override;
};


#endif  // SWITCH_ITERATION_CALLBACK
