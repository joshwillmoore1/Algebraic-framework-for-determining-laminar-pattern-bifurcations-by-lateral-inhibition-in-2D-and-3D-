#include "CollierSrnModel.hpp"

CollierSrnModel::CollierSrnModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractOdeSrnModel(2, pOdeSolver)
{
    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<CollierSrnModel, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        mpOdeSolver->SetMaxSteps(10000);
#else
        mpOdeSolver = CellCycleModelOdeSolver<MamDeltaNotchSrnModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
        SetDt(0.001);
#endif //CHASTE_CVODE
    }
    assert(mpOdeSolver->IsSetUp());
}

CollierSrnModel::CollierSrnModel(const CollierSrnModel& rModel)
    : AbstractOdeSrnModel(rModel)
{
    /*
     * Set each member variable of the new SRN model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new SRN model's member variables
     * will already have been correctly initialized in its constructor.
     *
     * Note 2: one or more of the new SRN model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new SRN model.
     *
     * Note 3: Only set the variables defined in this class. Variables defined
     * in parent classes will be defined there.
     */

    assert(rModel.GetOdeSystem());
    SetOdeSystem(new CollierOdeSystem(rModel.GetOdeSystem()->rGetStateVariables()));
}

AbstractSrnModel* CollierSrnModel::CreateSrnModel()
{
    return new CollierSrnModel(*this);
}

void CollierSrnModel::SimulateToCurrentTime()
{
    // Custom behaviour
    UpdateCollier();

    // Run the ODE simulation as needed
    AbstractOdeSrnModel::SimulateToCurrentTime();
}

void CollierSrnModel::Initialise()
{
    AbstractOdeSrnModel::Initialise(new CollierOdeSystem);
}

void CollierSrnModel::UpdateCollier()
{
    assert(mpOdeSystem != nullptr);
    assert(mpCell != nullptr);

    double mean_delta = mpCell->GetCellData()->GetItem("mean delta");
    mpOdeSystem->SetParameter("Mean Delta", mean_delta);

}

double CollierSrnModel::GetNotch()
{
    assert(mpOdeSystem != nullptr);
    double notch = mpOdeSystem->rGetStateVariables()[0];
    return notch;
}

double CollierSrnModel::GetDelta()
{
    assert(mpOdeSystem != nullptr);
    double delta = mpOdeSystem->rGetStateVariables()[1];
    return delta;
}

double CollierSrnModel::GetMeanNeighbouringDelta()
{
    assert(mpOdeSystem != nullptr);
    double mean_neighbouring_delta = mpOdeSystem->GetParameter("Mean Delta");
    return mean_neighbouring_delta;
}

void CollierSrnModel::OutputSrnModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractOdeSrnModel::OutputSrnModelParameters(rParamsFile);
}

// Declare identifier for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CollierSrnModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(CollierSrnModel)
