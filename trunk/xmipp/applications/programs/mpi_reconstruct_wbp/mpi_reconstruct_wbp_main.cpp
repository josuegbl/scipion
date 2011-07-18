/***************************************************************************
 *
 * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include <parallel/xmipp_mpi.h>
#include <reconstruction/reconstruct_wbp.h>

class ProgMPIRecWbp: public ProgRecWbp, public MpiMetadataProgram
{
public:
    void defineParams()
    {
        ProgRecWbp::defineParams();
        MpiMetadataProgram::defineParams();
    }
    void readParams()
    {
        MpiMetadataProgram::readParams();
        ProgRecWbp::readParams();
    }
    void read(int argc, char **argv)
    {
        MpiMetadataProgram::read(argc,argv);
        ProgRecWbp::read(argc, argv);
    }
    void produceSideInfo()
    {
        ProgRecWbp::produceSideInfo();
        createTaskDistributor(SF);
    }
    void showProgress()
    {
        if (node->isMaster() && verbose>0)
        {
            time_bar_done=first+1;
            progress_bar(time_bar_done);
        }
    }
    bool getImageToProcess(size_t &objId, size_t &objIndex)
    {
        return getTaskToProcess(objId, objIndex);
    }
    void finishProcessing()
    {
        MultidimArray<double> aux;
        aux.resizeNoCopy(reconstructedVolume());
        MPI_Allreduce(MULTIDIM_ARRAY(reconstructedVolume()), MULTIDIM_ARRAY(aux),
                      MULTIDIM_SIZE(aux), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        int iaux;
        MPI_Allreduce(&count_thr, &iaux, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        if (node->isMaster())
        {
            reconstructedVolume()=aux;
            count_thr=iaux;
            ProgRecWbp::finishProcessing();
        }
    }
};

int main(int argc, char **argv)
{
    ProgMPIRecWbp program;
    program.read(argc,argv);
    return program.tryRun();
}
