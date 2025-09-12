#include <TRandom3.h>
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <iostream>

int test() {
    // 定义随机数生成器
    TRandom3 randomGen(0); // 使用时间种子

    // 定义高斯分布的PDF
    TF1 gaussPDF("gaussPDF", "gaus", -5, 5);
    gaussPDF.SetParameters(1, 0, 1); // 参数：振幅=1，均值=0，标准差=1

    // 生成高斯的toy histogram
    int nBins = 50;
    TH1F toyHist("toyHist", "Gaussian Toy Histogram", nBins, -1, 1);
    for (int i = 0; i < 100; ++i) {
        toyHist.Fill(randomGen.Gaus(0, 1.0)); // 均值=0，标准差=1
    }

    // 生成PDF的histogram
    TH1F pdfHist("pdfHist", "PDF Histogram", nBins, -1, 1);
    for (int i = 1; i <= nBins; ++i) {
        double binCenter = pdfHist.GetBinCenter(i);
        double binContent = gaussPDF.Eval(binCenter);
        pdfHist.SetBinContent(i, binContent);
    }

    // 归一化两个直方图
    toyHist.Scale(1.0 / toyHist.Integral());
    pdfHist.Scale(1.0 / pdfHist.Integral());

    // 使用KolmogorovTest计算prob
    double prob = toyHist.KolmogorovTest(&pdfHist); // Use the TH1::KolmogorovTest method
    std::cout << "Kolmogorov Test Probability: " << prob << std::endl;

    return 0;
}