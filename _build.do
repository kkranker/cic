*! _build.do
*  Prep files for release

version 15.1
translator set smcl2txt linesize 120

translate "C:\Users\kkranker\OneDrive - Mathematica\Documents\Stata\cic\code-cic\cic.sthlp"  ///
          "C:\Users\kkranker\OneDrive - Mathematica\Documents\Stata\cic\code-cic\cic.txt"   ///
          , replace translator(smcl2txt)
