#include "pdfXY.h"

namespace MTMCSim {


//should not be too large like 1e-5, otw will ignore some terms of joint entropy


void pdfXY::compEntro()
{

    HX = 0;
    double ln2 = log(2.0);
    for (unsigned int i = 0; i<abXSize; i++)
        if (xMarPdf[i] <= EPS)
            continue;
        else
        HX +=  xMarPdf[i] * log(1/xMarPdf[i]) / ln2;

    HY = 0;
    for (unsigned int i = 0; i<abYSize; i++)
        if (yMarPdf[i] <= EPS)
            continue;
        else

        HY += yMarPdf[i] * log(1/yMarPdf[i]) / ln2;

    joH = 0;
    for (unsigned int i = 0; i<abXSize; i++)
        for (unsigned int j = 0; j<abYSize; j++)
            if (joPdf[i][j] <= EPS)
                continue;
            else
            joH += joPdf[i][j] * log(1/joPdf[i][j]) / ln2;

    condHYX = joH - HX;
    condHXY = joH - HY;
    IXY = HX + HY - joH;
}

void pdfXY::setPdf(const std::vector<double>& xMarPdf_, const std::vector<vector<double> >& yxCondPdf_)
{
    abXSize = xMarPdf_.size();
    abYSize = yxCondPdf_[0].size(); // assume each row has same # of cols
    

    joPdf.resize(abXSize);

    for (unsigned int i = 0; i<abXSize; i++)
        joPdf[i].resize(abYSize);

    xMarPdf = xMarPdf_;
    yxCondPdf = yxCondPdf_;

    joPdf.resize(abXSize);

    for (unsigned int i = 0; i<abXSize; i++)
        joPdf[i].resize(abYSize);


    for (unsigned int i = 0; i<abXSize; i++)
        for (unsigned int j = 0; j<abYSize; j++)
            joPdf[i][j] = xMarPdf[i] * yxCondPdf[i][j];


    yMarPdf.resize(abYSize);

    for (unsigned int i = 0; i<abYSize; i++)
        yMarPdf[i] = 0;

    for (unsigned int i = 0; i<abXSize; i++)
        for (unsigned int j = 0; j<abYSize; j++)
            yMarPdf[j] += joPdf[i][j];

    compEntro();
}

void pdfXY::setPdf(const std::vector<vector<double> > &joPdf_)
{
    abXSize = joPdf_.size();
    abYSize = joPdf_[0].size();


    joPdf = joPdf_;

    xMarPdf.resize(abXSize);
    yMarPdf.resize(abYSize);
    yxCondPdf.resize(abXSize);
    for (unsigned int i = 0; i<abXSize; i++)
        yxCondPdf[i].resize(abYSize);

    for (unsigned int i = 0; i<abXSize; i++)
        xMarPdf[i] = 0;

    for (unsigned int i = 0; i<abYSize; i++)
        yMarPdf[i] = 0;

    for (unsigned int i = 0; i<abXSize; i++)
        for (unsigned int j = 0; j<abYSize; j++)
        {
            xMarPdf[i] += joPdf[i][j];
            yMarPdf[j] += joPdf[i][j];
        }

    for (unsigned int i = 0; i<abXSize; i++)
        for (unsigned int j = 0; j<abYSize; j++)
        {
            if (xMarPdf[i] >= EPS)
                yxCondPdf[i][j] = (double)joPdf[i][j] / xMarPdf[i];
            else
                yxCondPdf[i][j] = 0;
        }

    compEntro();
}



void pdfXY::setPdf(double agreeProb, int size)
{
    abXSize = size;
    abYSize = size;

    xMarPdf.resize(abXSize);

    for (unsigned int i = 0; i<abXSize; i++)
        xMarPdf[i] = (double)1/abXSize;

    yMarPdf.resize(abYSize);

    for (unsigned int i = 0; i<abYSize; i++)
        yMarPdf[i] = (double)1/abYSize;

    yxCondPdf.resize(abXSize);

    for (unsigned int i = 0; i<abXSize; i++)
    {
        yxCondPdf[i].resize(abYSize);
        for (unsigned int j = 0; j<abYSize; j++)
        {
            if (i == j)
                yxCondPdf[i][j] = agreeProb;
            else
                yxCondPdf[i][j] = (double)(1-agreeProb)/(abYSize - 1);
        }
    }

    joPdf.resize(abXSize);

    for (unsigned int i = 0; i<abXSize; i++)
    {
        joPdf[i].resize(abYSize);
        
        for (unsigned int j = 0; j<abYSize; j++)
            joPdf[i][j] = xMarPdf[i] * yxCondPdf[i][j];
    }

    compEntro();
}

void pdfXY::setPdf(std::string pdffn, int xSize, int ySize)
{
    abXSize = xSize;
    abYSize = ySize;

    ifstream pdfFile(pdffn.c_str());

    if (!pdfFile)
        throw runtime_error("Read File Error!" + pdffn);

    joPdf.resize(abXSize);

    for (unsigned int i = 0; i<abXSize; i++)
    {
        joPdf[i].resize(abYSize);
        for (unsigned int j = 0; j<abYSize; j++)
            if (!(pdfFile >> joPdf[i][j]))
                throw runtime_error("PDF File Input Error!");
    }

    pdfFile.close();

    xMarPdf.resize(abXSize);
    yMarPdf.resize(abYSize);
    yxCondPdf.resize(abXSize);
    for (unsigned int i = 0; i<abXSize; i++)
        yxCondPdf[i].resize(abYSize);

    for (unsigned int i = 0; i<abXSize; i++)
        xMarPdf[i] = 0;

    for (unsigned int i = 0; i<abYSize; i++)
        yMarPdf[i] = 0;

    for (unsigned int i = 0; i<abXSize; i++)
        for (unsigned int j = 0; j<abYSize; j++)
        {
            xMarPdf[i] += joPdf[i][j];
            yMarPdf[j] += joPdf[i][j];
        }

    for (unsigned int i = 0; i<abXSize; i++)
        for (unsigned int j = 0; j<abYSize; j++)
        {
            if (xMarPdf[i] >= EPS)
                yxCondPdf[i][j] = (double)joPdf[i][j] / xMarPdf[i];
            else
                yxCondPdf[i][j] = 0;
        }

    compEntro();
}








}

//int main()
//{
//  using namespace MTMCSim;
//  pdfXY qary09_256(0.9, 256);
//  pdfXY ffile("Test/SPQ256NZ8Ran_1.txt", 256,256);
//}
