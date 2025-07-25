/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________    ______ ______ _    _
                       | | | ||  ___|  _  \_   _| ___ \   |  _  \|  ___| \  / |
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /   | | | || |_  |  \/  |
 / _ \| '_ \ / _ \ '_ \|  _  ||  _| | | | | | | | ___ \---| | | ||  _| | |\/| |
| (_) | |_) |  __/ | | | | | || |   | |/ / _| |_| |_/ /---| |/ / | |___| |  | |
 \___/| .__/ \___|_| |_\_| |_/\_|   |___/  \___/\____/    |___/  |_____|_|  |_|
      | |                     H ybrid F ictitious D omain - I mmersed B oundary
      |_|                                        and D iscrete E lement M ethod
-------------------------------------------------------------------------------
License

    openHFDIB-DEM is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with openHFDIB. If not, see <http://www.gnu.org/licenses/lgpl.html>.

InNamspace
    Foam

Description
    class for interpolation

SourceFiles

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef interpolationInfo_H
#define interpolationInfo_H

#include "dictionary.H"
#include "fvCFD.H"
// #include "fvMesh.H"
#include "triSurfaceMesh.H"
#include "geomModel.H"

namespace Foam
{
    class interpolationInfo
    {
        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

        protected:

        // mesh reference
        const Foam::fvMesh& mesh_;
        // geomModel reference
        std::shared_ptr<geomModel>& geomModel_;

        autoPtr<List<point>> ibPointsPtr_; // cellIbPoint mapping
        autoPtr<List<vector>> ibNormalsPtr_; // cellIbNormal mapping

        public:

        // Constructors ------------------------------------------------
        interpolationInfo
        (
            const Foam::fvMesh& mesh,
            std::shared_ptr<geomModel>& gModel
        );

        // Destructors -------------------------------------------------
        virtual ~interpolationInfo();

        List<point>& getIbPoints()
        {
            if(ibPointsPtr_->empty())
            {
                setIntpInfo();
            }
            return ibPointsPtr_();
        }

        List<vector>& getIbNormals()
        {
            return ibNormalsPtr_();
        }

        virtual void setIntpInfo() = 0;

        const DynamicLabelList& getSurfCells() const
        {
            return geomModel_->getSurfaceCellList()[Pstream::myProcNo()];
        }

        virtual void clearIntpInfo()
        {
            ibPointsPtr_.clear();
            ibNormalsPtr_.clear();
        }
    };
}
#endif
