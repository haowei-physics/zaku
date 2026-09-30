#include <TH2D.h>

#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include "zaku/stages/Stages.h"

namespace zaku {

RateTable::RateTable(std::vector<std::string> rows, std::vector<std::string> columns)
    : rows_(std::move(rows)), columns_(std::move(columns)) {
  weights_.assign(rows_.size() * columns_.size(), 0.0);
  counts_.assign(rows_.size() * columns_.size(), 0);
}

void RateTable::add(int row, int column, double weight) {
  if (row < 0 || row >= static_cast<int>(rows_.size()) || column < 0 || column >= static_cast<int>(columns_.size())) {
    throw std::out_of_range("RateTable::add: index out of range");
  }
  weights_[at(row, column)] += weight;
  ++counts_[at(row, column)];
}

double RateTable::weight(int row, int column) const { return weights_[at(row, column)]; }
long long RateTable::entries(int row, int column) const { return counts_[at(row, column)]; }

double RateTable::rowTotal(int row) const {
  double s = 0.0;
  for (std::size_t c = 0; c < columns_.size(); ++c) s += weights_[at(row, static_cast<int>(c))];
  return s;
}

long long RateTable::rowEntries(int row) const {
  long long n = 0;
  for (std::size_t c = 0; c < columns_.size(); ++c) n += counts_[at(row, static_cast<int>(c))];
  return n;
}

double RateTable::columnTotal(int column) const {
  double s = 0.0;
  for (std::size_t r = 0; r < rows_.size(); ++r) s += weights_[at(static_cast<int>(r), column)];
  return s;
}

double RateTable::total() const {
  double s = 0.0;
  for (double w : weights_) s += w;
  return s;
}

std::string RateTable::csv() const {
  std::ostringstream os;
  os << std::setprecision(10);
  os << "sample";
  for (const auto& c : columns_) os << "," << c;
  os << ",total,entries\n";
  for (std::size_t r = 0; r < rows_.size(); ++r) {
    os << rows_[r];
    for (std::size_t c = 0; c < columns_.size(); ++c) os << "," << weights_[at(static_cast<int>(r), static_cast<int>(c))];
    os << "," << rowTotal(static_cast<int>(r)) << "," << rowEntries(static_cast<int>(r)) << "\n";
  }
  return os.str();
}

std::string RateTable::text() const {
  std::ostringstream os;
  std::size_t w = 12;
  for (const auto& r : rows_) w = std::max(w, r.size() + 2);
  os << std::left << std::setw(static_cast<int>(w)) << "";
  for (const auto& c : columns_) os << std::right << std::setw(static_cast<int>(std::max<std::size_t>(14, c.size() + 2))) << c;
  os << std::setw(14) << "total" << std::setw(12) << "entries" << "\n";
  for (std::size_t r = 0; r < rows_.size(); ++r) {
    os << std::left << std::setw(static_cast<int>(w)) << rows_[r];
    for (std::size_t c = 0; c < columns_.size(); ++c) {
      os << std::right << std::setw(static_cast<int>(std::max<std::size_t>(14, columns_[c].size() + 2))) << std::fixed
         << std::setprecision(1) << weights_[at(static_cast<int>(r), static_cast<int>(c))];
    }
    os << std::setw(14) << std::fixed << std::setprecision(1) << rowTotal(static_cast<int>(r)) << std::setw(12)
       << rowEntries(static_cast<int>(r)) << "\n";
  }
  return os.str();
}

TH2D* RateTable::histogram(const std::string& name) const {
  auto* h = new TH2D(name.c_str(), "expected events;production;sample", static_cast<int>(columns_.size()), 0.0,
                     static_cast<double>(columns_.size()), static_cast<int>(rows_.size()), 0.0,
                     static_cast<double>(rows_.size()));
  h->SetDirectory(nullptr);
  for (std::size_t c = 0; c < columns_.size(); ++c) h->GetXaxis()->SetBinLabel(static_cast<int>(c) + 1, columns_[c].c_str());
  for (std::size_t r = 0; r < rows_.size(); ++r) {
    h->GetYaxis()->SetBinLabel(static_cast<int>(r) + 1, rows_[r].c_str());
    for (std::size_t c = 0; c < columns_.size(); ++c) {
      h->SetBinContent(static_cast<int>(c) + 1, static_cast<int>(r) + 1, weights_[at(static_cast<int>(r), static_cast<int>(c))]);
    }
  }
  return h;
}

RateTable RateTable::fromCsv(const std::string& text) {
  std::istringstream in(text);
  std::string line;
  if (!std::getline(in, line)) throw std::runtime_error("RateTable::fromCsv: empty text");
  std::vector<std::string> header;
  {
    std::istringstream ls(line);
    std::string cell;
    while (std::getline(ls, cell, ',')) header.push_back(cell);
  }
  if (header.size() < 3 || header.front() != "sample") throw std::runtime_error("RateTable::fromCsv: bad header");
  std::vector<std::string> columns(header.begin() + 1, header.end() - 2);
  std::vector<std::string> rows;
  std::vector<std::vector<double>> w;
  std::vector<long long> n;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    std::istringstream ls(line);
    std::string cell;
    std::vector<std::string> cells;
    while (std::getline(ls, cell, ',')) cells.push_back(cell);
    if (cells.size() != header.size()) throw std::runtime_error("RateTable::fromCsv: bad row '" + line + "'");
    rows.push_back(cells[0]);
    std::vector<double> ws;
    for (std::size_t c = 1; c + 2 < cells.size(); ++c) ws.push_back(std::stod(cells[c]));
    w.push_back(ws);
    n.push_back(std::stoll(cells.back()));
  }
  RateTable t(rows, columns);
  for (std::size_t r = 0; r < rows.size(); ++r) {
    for (std::size_t c = 0; c < columns.size(); ++c) t.weights_[t.at(static_cast<int>(r), static_cast<int>(c))] = w[r][c];
    if (!columns.empty()) t.counts_[t.at(static_cast<int>(r), 0)] = n[r];  // entries are only known per row
  }
  return t;
}

void RateTable::writeCsv(const std::string& path) const {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("RateTable::writeCsv: cannot write '" + path + "'");
  out << csv();
}

}  // namespace zaku
