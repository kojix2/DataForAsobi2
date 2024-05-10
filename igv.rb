require "htslib"
require "igv"

igv = IGV.start
igv.genome "ref.fa"
igv.load "s123_variants.vcf"
igv.load "s1.sorted.bam"
igv.load "s2.sorted.bam"
igv.load "s3.sorted.bam"
igv.squish

v = HTS::Bcf.open "s123_variants.vcf"
v.each do |r|
  puts r.to_s
  igv.go "#{r.chrom}:#{r.pos}"
  igv.snapshot "snapshot/#{r.chrom}_#{r.pos}.png"
  gets
end

v.close
igv.exit
