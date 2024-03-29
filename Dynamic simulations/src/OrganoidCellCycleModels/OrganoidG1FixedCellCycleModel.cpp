#include "OrganoidG1FixedCellCycleModel.hpp"
#include "BasalStemCellProliferativeType.hpp"
#include "MyoEpiCellProliferativeType.hpp"
#include "LumenERPositiveCellProliferativeType.hpp"
#include "LumenERNegativeCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"


OrganoidG1FixedCellCycleModel::OrganoidG1FixedCellCycleModel()
    : AbstractSimplePhaseBasedCellCycleModel()
{
}



OrganoidG1FixedCellCycleModel::OrganoidG1FixedCellCycleModel(const OrganoidG1FixedCellCycleModel& rModel)
   : AbstractSimplePhaseBasedCellCycleModel(rModel)
{
    /*
     * The member variables mGeneration and mMaxTransitGeneration are
     * initialized in the AbstractSimpleGenerationalCellCycleModel
     * constructor.
     *
     * The member variables mCurrentCellCyclePhase, mG1Duration,
     * mMinimumGapDuration, mStemCellG1Duration, mTransitCellG1Duration,
     * mSDuration, mG2Duration and mMDuration are initialized in the
     * AbstractPhaseBasedCellCycleModel constructor.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     *
     * Note that mG1Duration is (re)set as soon as InitialiseDaughterCell()
     * is called on the new cell-cycle model.
     * 
     * NOTE: THIS CLASS CURRENTLY HAS ARTIFICIAL CELL CYCLES TIMES
     * 
     */
}

AbstractCellCycleModel* OrganoidG1FixedCellCycleModel::CreateCellCycleModel()
{
    return new OrganoidG1FixedCellCycleModel(*this);
}

void OrganoidG1FixedCellCycleModel::SetG1Duration()
{

    //Stochastic G1 has been added -  using similar values as the Uniform generational cell cycle model
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    assert(mpCell != nullptr);

    //organoid cells are set to be diff for SRN investigation

    if (mpCell->GetCellProliferativeType()->IsType<BasalStemCellProliferativeType>())
    {
        mG1Duration = 18 + 16*p_gen->ranf(); 
    }
    else if (mpCell->GetCellProliferativeType()->IsType<MyoEpiCellProliferativeType>())
    {
        mG1Duration =  DBL_MAX;
    }
    else if (mpCell->GetCellProliferativeType()->IsType<LumenERPositiveCellProliferativeType>())
    {
        mG1Duration = DBL_MAX;
    }

    else if (mpCell->GetCellProliferativeType()->IsType<LumenERNegativeCellProliferativeType>())
    {
        mG1Duration = 6 + 8*p_gen->ranf();
    }

    else if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        mG1Duration = DBL_MAX;
    }

    else
    {   
        NEVER_REACHED;
    }
}

void OrganoidG1FixedCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractSimplePhaseBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(OrganoidG1FixedCellCycleModel)