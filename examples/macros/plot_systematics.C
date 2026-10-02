// Draw the systematics summary: for one sample, the nominal spectrum and the
// +-1 sigma spectra of the largest-effect parameters, plus a bar chart of the
// fractional effects.
//
//   root -l 'examples/macros/plot_systematics.C("work/summary/systematics_summary.root",
//                                               "work/summary/systematics_summary.csv", "FHC_1Rmu", 8)'
//
// Needs libZaku on LD_LIBRARY_PATH (source setup.sh) for zaku::SummaryTable.
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "zaku/Zaku.h"

void plot_systematics(const char* rootFile = "systematics_summary.root", const char* csvFile = "systematics_summary.csv",
                      const char* sample = "FHC_1Rmu", int nTop = 8) {
  gStyle->SetOptStat(0);
  const zaku::SummaryTable table = zaku::SummaryTable::readCsv(csvFile);

  // the largest effects on this sample
  std::vector<zaku::SummaryEntry> top;
  for (const auto& e : table.top(-1)) {
    if (e.sample == sample) top.push_back(e);
    if (static_cast<int>(top.size()) >= nTop) break;
  }

  TFile* f = TFile::Open(rootFile, "READ");
  TDirectory* dir = f->GetDirectory(sample);
  if (!dir) {
    printf("no sample '%s' in %s\n", sample, rootFile);
    return;
  }
  auto* c = new TCanvas("c_syst", sample, 1200, 500);
  c->Divide(2, 1);

  // left: spectra
  c->cd(1);
  auto* nominal = dir->Get<TH1D>("nominal");
  nominal->SetLineColor(kBlack);
  nominal->SetLineWidth(2);
  nominal->SetTitle((std::string(sample) + ";" + nominal->GetXaxis()->GetTitle() + ";events").c_str());
  nominal->Draw("hist");
  auto* leg = new TLegend(0.55, 0.55, 0.89, 0.89);
  leg->AddEntry(nominal, "nominal", "l");
  const int colours[] = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1};
  for (std::size_t i = 0; i < top.size() && i < 4; ++i) {
    auto* up = dir->Get<TH1D>((top[i].parameter + "_plus").c_str());
    auto* down = dir->Get<TH1D>((top[i].parameter + "_minus").c_str());
    if (!up || !down) continue;
    up->SetLineColor(colours[i]);
    down->SetLineColor(colours[i]);
    down->SetLineStyle(2);
    up->Draw("hist same");
    down->Draw("hist same");
    leg->AddEntry(up, (top[i].parameter + " #pm1#sigma").c_str(), "l");
  }
  leg->Draw();

  // right: fractional effects
  c->cd(2);
  auto* bars = new TH1D("bars", (std::string(sample) + " largest #pm1#sigma effects;;#Delta N / N").c_str(),
                        static_cast<int>(top.size()), 0, static_cast<double>(top.size()));
  auto* barsMinus = static_cast<TH1D*>(bars->Clone("bars_minus"));
  for (std::size_t i = 0; i < top.size(); ++i) {
    bars->GetXaxis()->SetBinLabel(static_cast<int>(i) + 1, top[i].parameter.c_str());
    bars->SetBinContent(static_cast<int>(i) + 1, top[i].plusFrac());
    barsMinus->SetBinContent(static_cast<int>(i) + 1, top[i].minusFrac());
  }
  bars->SetFillColor(kRed - 7);
  barsMinus->SetFillColor(kBlue - 7);
  const double range = 1.2 * std::max(std::abs(bars->GetMaximum()), std::abs(barsMinus->GetMinimum()));
  bars->SetMinimum(-range);
  bars->SetMaximum(range);
  bars->Draw("bar");
  barsMinus->Draw("bar same");
  bars->GetXaxis()->LabelsOption("v");
  c->Update();
}
