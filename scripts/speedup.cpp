// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <sstream>
using namespace Rcpp;

// ---------- helpers ----------
static inline std::vector<std::string> split(const std::string& s, char delim) {
  std::vector<std::string> out;
  std::string cur;
  std::istringstream iss(s);
  while (std::getline(iss, cur, delim)) out.push_back(cur);
  return out;
}

static inline std::string paste_c(const std::vector<std::string>& v, const std::string& sep) {
  if (v.empty()) return "";
  std::ostringstream oss;
  for (size_t i=0; i<v.size(); ++i) {
    if (i) oss << sep;
    oss << v[i];
  }
  return oss.str();
}

static inline bool parse_double_strict(const std::string& s, double& out) {
  char* end = nullptr;
  out = std::strtod(s.c_str(), &end);
  return end && *end == '\0';
}



// ---------- testContinue ----------
/*
  R: testContinue(vec) -> logical
  returns TRUE if numeric(vec) is strictly contiguous by +/-1 steps.
*/
// [[Rcpp::export]]
bool testContinueCpp(CharacterVector vec) {
  int n = vec.size();
  if (n <= 1) return true;
  std::vector<double> x(n);
  for (int i = 0; i < n; ++i) {
    if (vec[i] == NA_STRING) return false;
    std::string s = Rcpp::as<std::string>(vec[i]);
    double z;
    if (!parse_double_strict(s, z)) return false;
    x[i] = z;
  }
  for (int i = 1; i < n; ++i) {
    if (std::abs(x[i] - x[i - 1]) != 1.0) return false;
  }
  return true;
}

// ---------- isoformFilter ----------
/*
  R: isoformFilter(isoformData, tol = 9) -> logical
  isoformData columns used: start, end, CDS_start, CDS_end, strand
  NOTE: assumes column *start* (not 'strat') in allCDS; your original appears to have a typo.
*/
// [[Rcpp::export]]
bool isoformFilterCpp(DataFrame isoformData, int tol = 9) {
  IntegerVector start = isoformData["start"];
  IntegerVector end   = isoformData["end"];
  IntegerVector CDS_start = isoformData["CDS_start"];
  IntegerVector CDS_end   = isoformData["CDS_end"];
  CharacterVector strandV = isoformData["strand"];
  if (strandV.size() == 0) return false;

  std::string strand = Rcpp::as<std::string>(strandV[0]);
  int n = isoformData.nrows();

  // If strand is not '+' or '-', we cannot match edges reliably
  if (strand != "+" && strand != "-") return false;

  // Copy so we can reorder for '-' strand
  std::vector<int> s(start.begin(), start.end());
  std::vector<int> e(end.begin(), end.end());
  std::vector<int> cs(CDS_start.begin(), CDS_start.end());
  std::vector<int> ce(CDS_end.begin(), CDS_end.end());

  if (strand == "-") {
    std::reverse(s.begin(), s.end());
    std::reverse(e.begin(), e.end());
    std::reverse(cs.begin(), cs.end());
    std::reverse(ce.begin(), ce.end());
  }

  if (n == 1) {
    return (s[0] + tol >= cs[0]) && (e[0] - tol <= ce[0]);
  }

  if (strand == "+") {
    if (!((s[0] + tol >= cs[0]) && (std::abs(e[0] - ce[0]) <= tol))) return false;
    if (!((std::abs(s[n-1] - cs[n-1]) <= tol) && (e[n-1] - tol <= ce[n-1]))) return false;
    for (int i = 1; i < n-1; ++i)
      if (std::abs(s[i] - cs[i]) > tol || std::abs(e[i] - ce[i]) > tol) return false;
    return true;
  } else { // '-'
    if (!((e[0] - tol <= ce[0]) && (std::abs(s[0] - cs[0]) <= tol))) return false;
    if (!((std::abs(e[n-1] - ce[n-1]) <= tol) && (s[n-1] + tol >= cs[n-1]))) return false;
    for (int i = 1; i < n-1; ++i)
      if (std::abs(s[i] - cs[i]) > tol || std::abs(e[i] - ce[i]) > tol) return false;
    return true;
  }
}

// ---------- perGene ----------
/*
  perGene(geneData, allCDS, tol=9) returns a character vector of 8 fields:
  1 qual.source (||-joined)
  2 position (block keys for first isoform)
  3 nblock (as string)
  4 nCDS_ref (||-joined)
  5 length_ref (||-joined)
  6 note
  7 NO.CDS (||-joined order values)
  8 type

  Columns expected in geneData:
    chr, start, end, strand, SampleID, CDS_chr, CDS_start, CDS_end, CDS_strand,
    n_base, gene, isoform, order
  Columns expected in allCDS:
    isoform, start, end
*/
// [[Rcpp::export]]
CharacterVector perGeneCpp(DataFrame geneData, DataFrame allCDS, int tol = 9) {
  CharacterVector chr    = geneData["chr"];
  IntegerVector   start  = geneData["start"];
  IntegerVector   end    = geneData["end"];
  CharacterVector strand = geneData["strand"];
  CharacterVector iso    = geneData["isoform"];
  CharacterVector orderV = geneData["order"];
  IntegerVector   cdsS   = geneData["CDS_start"];
  IntegerVector   cdsE   = geneData["CDS_end"];

  const int n = geneData.nrows();

  // Build block keys and isoform -> row indices
  std::vector<std::string> blockKey(n);
  std::unordered_map<std::string, std::vector<int>> isoRows;
  blockKey.reserve(n);
  for (int i = 0; i < n; ++i) {
    std::string k = Rcpp::as<std::string>(chr[i]) + ":" +
                    std::to_string(start[i]) + ":" +
                    std::to_string(end[i]) + ":" +
                    Rcpp::as<std::string>(strand[i]);
    blockKey[i] = k;
    isoRows[Rcpp::as<std::string>(iso[i])].push_back(i);
  }

  // Unique blocks in this read
  std::vector<std::string> uniqBlocks = blockKey;
  std::sort(uniqBlocks.begin(), uniqBlocks.end());
  uniqBlocks.erase(std::unique(uniqBlocks.begin(), uniqBlocks.end()), uniqBlocks.end());
  const int nblock = (int)uniqBlocks.size();

  // Helper: rows for a given iso
  auto get_rows = [&](const std::string& isof)->std::vector<int>{
    auto it = isoRows.find(isof);
    if (it == isoRows.end()) return {};
    return it->second;
  };

  // Candidate isoforms that cover ALL blocks and have a unique strand
  std::vector<std::string> candidates;
  candidates.reserve(isoRows.size());
  for (auto &kv : isoRows) {
    const std::string &isof = kv.first;
    // Build set of blocks for this iso
    std::unordered_set<std::string> seen;
    for (int idx : kv.second) seen.insert(blockKey[idx]);
    bool covers_all = true;
    for (auto &b : uniqBlocks) if (!seen.count(b)) { covers_all = false; break; }
    if (!covers_all) continue;

    // Unique strand?
    std::unordered_set<std::string> sset;
    for (int idx : kv.second) sset.insert(Rcpp::as<std::string>(strand[idx]));
    if (sset.size() == 1) candidates.push_back(isof);
  }

  // "position" should follow the first candidate isoform row order (as in R),
  // falling back to all unique blocks if no candidate exists.
  auto build_position = [&](const std::string& isof)->std::string{
    auto rows = get_rows(isof);
    std::vector<std::string> v;
    v.reserve(rows.size());
    for (int idx : rows) {
      v.push_back(Rcpp::as<std::string>(chr[idx]) + ":" +
                  std::to_string(start[idx]) + ":" +
                  std::to_string(end[idx]) + ":" +
                  Rcpp::as<std::string>(strand[idx]));
    }
    return paste_c(v, ";");
  };
  std::string position_joined;
  if (!candidates.empty()) {
    position_joined = build_position(candidates[0]);
  } else {
    position_joined = paste_c(uniqBlocks, ";");
  }

  // Summaries from allCDS
  CharacterVector acIso = allCDS["isoform"];
  IntegerVector   acS   = allCDS["start"];   // NOTE: 'start' (fixes 'strat' typo)
  IntegerVector   acE   = allCDS["end"];

  auto ref_n_exon = [&](const std::string& isof)->int{
    int c = 0;
    for (int i=0;i<acIso.size();++i)
      if (acIso[i] != NA_STRING && Rcpp::as<std::string>(acIso[i]) == isof) ++c;
    return c;
  };
  auto ref_len = [&](const std::string& isof)->long long{
    long long L = 0;
    for (int i=0;i<acIso.size();++i)
      if (acIso[i] != NA_STRING && Rcpp::as<std::string>(acIso[i]) == isof)
        L += (long long)acE[i] - (long long)acS[i];
    return L;
  };

  // If there are candidates but all are "undefined" → novel with undefined
  if (!candidates.empty()) {
    bool all_undef = true;
    for (auto &c : candidates) if (c != "undefined") { all_undef = false; break; }
    if (all_undef) {
      std::vector<std::string> out(8);
      out[0] = paste_c(candidates, "||");           // all "undefined"
      out[1] = position_joined;
      out[2] = std::to_string(nblock);
      out[3] = "NA";
      out[4] = "NA";
      out[5] = "no full-covered isoform";
      out[6] = "NA";
      out[7] = "novel with undefined";
      return wrap(out);
    }
  }

  // No candidates at all → novel with addition
  if (candidates.empty()) {
    std::vector<std::string> out(8);
    out[0] = "undefined";
    out[1] = position_joined;
    out[2] = std::to_string(nblock);
    out[3] = "NA";
    out[4] = "NA";
    out[5] = "no full-covered isoform";
    out[6] = "NA";
    out[7] = "novel with addition";
    return wrap(out);
  }

  // Keep only continuous candidates (testContinue on 'order')
  auto orders_for_iso = [&](const std::string& isof)->CharacterVector{
    auto rows = get_rows(isof);
    CharacterVector v(rows.size());
    for (size_t i=0;i<rows.size();++i) v[i] = orderV[rows[i]];
    return v;
  };

  std::vector<std::string> cont_ok;
  for (const auto& isof : candidates) {
    if (testContinueCpp(orders_for_iso(isof))) cont_ok.push_back(isof);
  }

  if (cont_ok.empty()) {
    // discontinuous CDS → novel with deletion
    std::vector<std::string> NO_CDS, nCDS_s, len_s;
    for (const auto& isof : candidates) {
      auto rows = get_rows(isof);
      std::vector<std::string> ord; ord.reserve(rows.size());
      for (int idx : rows) ord.push_back(Rcpp::as<std::string>(orderV[idx]));
      NO_CDS.push_back(paste_c(ord, "-"));
      nCDS_s.push_back(std::to_string(ref_n_exon(isof)));
      len_s.push_back(std::to_string(ref_len(isof)));
    }
    std::vector<std::string> out(8);
    out[0] = paste_c(candidates, "||");
    out[1] = position_joined;
    out[2] = std::to_string(nblock);
    out[3] = paste_c(nCDS_s, "||");
    out[4] = paste_c(len_s, "||");
    out[5] = "discontinuous CDS";
    out[6] = paste_c(NO_CDS, "||");
    out[7] = "novel with deletion";
    return wrap(out);
  }

  // Edge matching on continuous ones
  auto df_for_iso = [&](const std::string& isof)->DataFrame{
    auto rows = get_rows(isof);
    const int m = (int)rows.size();
    IntegerVector s(m), e(m), cs(m), ce(m);
    CharacterVector st(m);
    for (int i=0;i<m;++i) {
      int idx = rows[i];
      s[i]  = start[idx];
      e[i]  = end[idx];
      cs[i] = cdsS[idx];
      ce[i] = cdsE[idx];
      st[i] = strand[idx];
    }
    return DataFrame::create(_["start"]=s,_["end"]=e,_["CDS_start"]=cs,_["CDS_end"]=ce,_["strand"]=st);
  };

  std::vector<std::string> edge_ok;
  for (const auto& isof : cont_ok) {
    if (isoformFilterCpp(df_for_iso(isof), tol)) edge_ok.push_back(isof);
  }

  const bool have_edge_ok = !edge_ok.empty();
  const std::vector<std::string>& chosen = have_edge_ok ? edge_ok : cont_ok;
  const std::string note = have_edge_ok ? "continuous CDS and edge-matching"
                                        : "continuous CDS but edge-unmatching";
  const std::string typ  = have_edge_ok ? "normal" : "novel with deletion";

  std::vector<std::string> NO_CDS, nCDS_s, len_s;
  for (const auto& isof : chosen) {
    auto rows = get_rows(isof);
    std::vector<std::string> ord; ord.reserve(rows.size());
    for (int idx : rows) ord.push_back(Rcpp::as<std::string>(orderV[idx]));
    NO_CDS.push_back(paste_c(ord, "-"));
    nCDS_s.push_back(std::to_string(ref_n_exon(isof)));
    len_s.push_back(std::to_string(ref_len(isof)));
  }

  std::vector<std::string> out(8);
  out[0] = paste_c(chosen, "||");
  out[1] = position_joined;
  out[2] = std::to_string(nblock);
  out[3] = paste_c(nCDS_s, "||");
  out[4] = paste_c(len_s, "||");
  out[5] = note;
  out[6] = paste_c(NO_CDS, "||");
  out[7] = typ;
  return wrap(out);
}

// ---------- filtGene ----------
/*
  filtGene(temp.info) -> list(outGene=list of character vectors, geneNum, nblock, position, overlapping)
  Re-implements the slow combn search with early-exit. For many genes we cap exhaustive search to size 8 and fall back to a greedy.
*/
// [[Rcpp::export]]
List filtGeneCpp(DataFrame info) {
  CharacterVector chr = info["chr"];
  IntegerVector start = info["start"];
  IntegerVector end   = info["end"];
  CharacterVector strand = info["strand"];
  CharacterVector gene = info["gene"];

  int n = info.nrows();
  std::vector<std::string> blocks(n);
  for (int i=0;i<n;++i) {
    blocks[i] = Rcpp::as<std::string>(chr[i]) + ":" + std::to_string(start[i]) + ":" + std::to_string(end[i]) + ":" + Rcpp::as<std::string>(strand[i]);
  }
  // Unique blocks + position string
  std::vector<std::string> uniqBlocks = blocks;
  std::sort(uniqBlocks.begin(), uniqBlocks.end());
  uniqBlocks.erase(std::unique(uniqBlocks.begin(), uniqBlocks.end()), uniqBlocks.end());
  std::string position = paste_c(uniqBlocks, ";");
  int nblock = (int)uniqBlocks.size();

  // Map gene -> set of covered blocks
  std::unordered_map<std::string, std::unordered_set<std::string>> G;
  for (int i=0;i<n;++i) {
    G[Rcpp::as<std::string>(gene[i])].insert(blocks[i]);
  }
  // Collect gene names
  std::vector<std::string> genes;
  genes.reserve(G.size());
  for (auto& kv : G) genes.push_back(kv.first);
  int m = (int)genes.size();

  // Build an index for blocks
  std::unordered_map<std::string,int> bIndex;
  for (int i=0;i<nblock;++i) bIndex[uniqBlocks[i]] = i;

  // Build binary coverage matrix (as vectors of bool)
  std::vector<std::vector<uint8_t>> cov(m, std::vector<uint8_t>(nblock, 0));
  for (int g=0; g<m; ++g) {
    for (auto& b : G[genes[g]]) {
      auto it = bIndex.find(b);
      if (it != bIndex.end()) cov[g][it->second] = 1;
    }
  }

  // Helper lambdas to test a combination
  auto covers_exactly_once = [&](const std::vector<int>& idxs)->bool{
    for (int b=0;b<nblock;++b) {
      int s=0; for (int g : idxs) s += cov[g][b];
      if (s != 1) return false;
    }
    return true;
  };
  auto covers_at_least_once = [&](const std::vector<int>& idxs)->bool{
    for (int b=0;b<nblock;++b) {
      int s=0; for (int g : idxs) s += cov[g][b];
      if (s < 1) return false;
    }
    return true;
  };

  // Exhaustive search with increasing size, up to 8 genes; then greedy fallback.
  auto search = [&](bool exact)->std::vector<std::vector<int>>{
    std::vector<std::vector<int>> res;
    int cap = std::min(m, 8); // cap exhaustive
    for (int k=1; k<=cap && res.empty(); ++k) {
      // iterate combinations of size k
      std::vector<int> idx(k);
      std::iota(idx.begin(), idx.end(), 0);
      while (true) {
        bool ok = exact ? covers_exactly_once(idx) : covers_at_least_once(idx);
        if (ok) res.push_back(idx);
        // next combination
        int pos = k - 1;
        while (pos >= 0 && idx[pos] == (m - k + pos)) --pos;
        if (pos < 0) break;
        ++idx[pos];
        for (int j=pos+1; j<k; ++j) idx[j] = idx[j-1] + 1;
      }
    }
    if (!res.empty() || m<=8) return res;

    // Greedy fallback (at least once or exactly once depending on 'exact')
    std::vector<int> chosen;
    std::vector<int> cover(nblock, 0);
    std::vector<int> all(m); std::iota(all.begin(), all.end(), 0);
    while ((int)chosen.size() < m) {
      int best = -1, bestGain = -1;
      for (int g : all) {
        if (std::find(chosen.begin(), chosen.end(), g) != chosen.end()) continue;
        int gain = 0;
        for (int b=0;b<nblock;++b) {
          int after = cover[b] + (cov[g][b] ? 1 : 0);
          if (exact) gain += (after == 1) && (cover[b] == 0);
          else       gain += (after >= 1) && (cover[b] == 0);
        }
        if (gain > bestGain) { bestGain = gain; best = g; }
      }
      if (best < 0) break;
      chosen.push_back(best);
      for (int b=0;b<nblock;++b) cover[b] += cov[best][b];
      // Early success check
      bool ok = true;
      for (int b=0;b<nblock;++b) {
        if (exact) { if (cover[b] != 1) { ok=false; break; } }
        else       { if (cover[b] < 1)  { ok=false; break; } }
      }
      if (ok) { res.push_back(chosen); break; }
    }
    return res;
  };

  auto exact = search(true);
  std::string overlapping = "N";
  std::vector<std::vector<int>> combos;
  if (!exact.empty()) combos = exact;
  else {
    overlapping = "Y";
    combos = search(false);
  }

  // Convert to list of gene name vectors
  List outGene;
  for (auto &c : combos) {
    CharacterVector v;
    for (int g : c) v.push_back(genes[g]);
    outGene.push_back(v);
  }
  return List::create(
    _["outGene"] = outGene,
    _["geneNum"] = combos.empty() ? 0 : (int)Rcpp::as<CharacterVector>(outGene[0]).size(),
    _["nblock"]  = nblock,
    _["position"]= position,
    _["overlapping"] = overlapping
  );
}
