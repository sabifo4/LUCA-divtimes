#! /usr/bin/env ruby


############################################################
require 'getoptlong'
require 'parallel'
require 'fileutils'
require 'tmpdir'
require 'find'
require 'colorize'

$LOAD_PATH << './lib'

require 'Dir.rb'
require 'processbar.rb'
require 'do_mcmctree.rb'


############################################################
#RSCRIPT="/usr/bin/Rscript" # otherwise cannot work w/ cl007
RSCRIPT='Rscript'

DIR = File.dirname(__FILE__)
LIB_DIR = File.join(DIR, 'lib')

IQTREE = 'iqtree2'
NW_STATS = 'nw_stats'
NW_TOPOLOGY = 'nw_topology'
MCMCTREE = 'mcmctree'
REORDER_NODE = File.join(DIR, 'reorder_node.rb')
FROM_BS_TO_HESSIAN = File.join(DIR, 'from_bs_to_hessian.R')
FIGTREE2NWK = File.expand_path(File.join(LIB_DIR, 'figtree2tree.sh'))


############################################################
def help()
  STDERR.puts "help message"
  exit 1
end


def make_argu(argu, value)
  return([argu, value].join(' '))
end


def processbar_for_bootstrapping(file:, b:)
  thr = Thread.new do |i|
    count = 0
    while true do
      next if not File.exist?(file)
      new_count = `wc -l #{file} | awk '{print $1}'`.chomp.to_i
      sleep(0.2)
      if new_count != count
        count = new_count
        processbar(count, b)
      end
    end
  end
  return(thr)
end


def get_files_under_folder(folder, b)
  files = Array.new
  Find.find(folder) do |path|
    files << path if path =~ /#{b}$/
  end
  return(files.join(' '))
end


def create_inBV(mltree_file:, mcmctree_outdir:, inBV_file:, iqtree_outdir:)
  no_species = `#{NW_STATS} #{mltree_file} | grep '^#leaves:' | awk '{print $2}'`.chomp.to_i  
  `bash -c "echo -e '\n#{no_species}\n' " >#{inBV_file}`

  `#{NW_TOPOLOGY} -bI #{mltree_file} >> #{inBV_file}`
  `echo >> #{inBV_file}`

  `cat #{iqtree_outdir}/ml.bls >> #{inBV_file}`
  `echo >> #{inBV_file}`

  gradient = [%w[0] * (2*no_species-3)].join(' ') # no. of branches equals 2n-3 where n is the no. of species
  `echo #{gradient} >> #{inBV_file}`
  `echo >> #{inBV_file}`

  `echo Hessian >> #{inBV_file}`
  `echo >> #{inBV_file}`
  `cat #{iqtree_outdir}/hessian >> #{inBV_file}`
end


def get_no_iter(mcmctree_ctl)
  no = nil
  in_fh = File.open(mcmctree_ctl, 'r')
  in_fh.each_line do |line|
    line.chomp!
    no = $1.to_i if line =~ /nsample[ ]*=[ ]*(\d+)/
  end
  in_fh.close
  return(no)
end


def split_ali_file(ali_file)
  in_fh = File.open(ali_file)
  ali2lines = Hash.new
  count = -1
  in_fh.each_line do |line|
    line.chomp!
    if line =~ /\d+\s+\d+$/
      count += 1
      ali2lines[count] = Array.new
    end
    ali2lines[count] << line
  end
  in_fh.close
  return(ali2lines)
end


############################################################
mcmctree_indir = nil
mcmctree_ctl_file = nil
ali_file = nil
model_argu = '-m LG+G'
bootstrap = 1000
bootstrap_argu = "-b #{bootstrap}"
te_argu = nil
is_pmsf = false
add_argu = nil
calib_tree_file = nil
ref_tree_file = nil

outdir = nil
is_force = false
cpu = 1
is_run_mcmctree = false


############################################################
opts = GetoptLong.new(
  ['--mcmctree_indir', GetoptLong::REQUIRED_ARGUMENT],
  ['--mcmctree_ctl', GetoptLong::REQUIRED_ARGUMENT],
  ['--ali', GetoptLong::REQUIRED_ARGUMENT],
  ['--ref', '--ref_tree_file', GetoptLong::REQUIRED_ARGUMENT],
  ['-m', GetoptLong::REQUIRED_ARGUMENT],
  ['-b', GetoptLong::REQUIRED_ARGUMENT],
  ['--pmsf', GetoptLong::NO_ARGUMENT],
  ['--te', GetoptLong::REQUIRED_ARGUMENT],
  ['--calib_tree', '--calibrated_tree', GetoptLong::REQUIRED_ARGUMENT],
  ['--outdir', GetoptLong::REQUIRED_ARGUMENT],
  ['--force', GetoptLong::NO_ARGUMENT],
  ['--cpu', GetoptLong::REQUIRED_ARGUMENT],
  ['--run_mcmctree', GetoptLong::NO_ARGUMENT],
  ['--add_cmd', '--add_argu', GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when '--mcmctree_indir'
      mcmctree_indir = value
    when '--mcmctree_ctl'
      mcmctree_ctl_file = value
    when '--ali'
      ali_file = value
    when '--ref', '--ref_tree_file'
      ref_tree_file = value
      te_argu = make_argu('-te', value)
    when '-m'
      model_argu = make_argu('-m', value)
    when '-b'
      bootstrap = value.to_i
      bootstrap_argu = make_argu('-b', value)
    when '--te'
      te_argu = make_argu('-te', value)
    when '--pmsf'
      is_pmsf = true
    when '--calib_tree', '--calibrated_tree'
      calib_tree_file = value
    when '--outdir'
      outdir = value
    when '--force'
      is_force = true
    when '--cpu'
      cpu = value.to_i
    when '--run_mcmctree'
      is_run_mcmctree = true
    when '-h'
      help()
    when '--add_cmd', '--add_argu'
      add_argu = value
  end
end


############################################################
mkdir_with_force(outdir, is_force)

ali2lines = split_ali_file(ali_file)

outdir_split = File.join(outdir, "split")

mkdir_with_force(outdir_split, is_force)

STDERR.puts "CPU:\t#{cpu}".colorize(:red)


############################################################
ali2lines.each_pair do |count, lines|
  outdir_split1 = File.join(outdir_split, count.to_s)
  mkdir_with_force(outdir_split1)
  ali_file1 = File.join(outdir_split1, 'combined.phy')
  out_fh = File.open(ali_file1, 'w')
  lines.each do |line|
    out_fh.puts line
  end
  out_fh.close

  iqtree_outdir = File.join(outdir_split1, 'iqtree')
  mcmctree_outdir = File.join(outdir_split1, 'mcmctree')

  mkdir_with_force(iqtree_outdir)
  mkdir_with_force(mcmctree_outdir)

  boottree_file = File.join(iqtree_outdir, 'iqtree.boottrees')
  mltree_file = File.join(iqtree_outdir, 'iqtree.treefile')
  inBV_file = File.join(mcmctree_outdir, 'in.BV')
  STDOUT.puts "Running IQ-Tree ......"

  # to record how many bootstrap trees are built
  thr = processbar_for_bootstrapping(file:boottree_file, b:bootstrap)

  if is_pmsf
    `#{IQTREE} -redo -s #{ali_file1} -pre #{iqtree_outdir}/guide -T #{cpu} -quiet -m LG4M+G #{te_argu} #{add_argu}`
    `#{IQTREE} -redo -s #{ali_file1} -pre #{iqtree_outdir}/iqtree -T #{cpu} -quiet #{model_argu} #{bootstrap_argu} #{te_argu} #{add_argu} -ft #{iqtree_outdir}/guide.treefile -blmin 1e-10`
  else
    STDERR.puts "pmsf not used".colorize(:blue) if model_argu =~ /C[0-9]+/
    `#{IQTREE} -redo -s #{ali_file1} -pre #{iqtree_outdir}/iqtree -T #{cpu} -quiet #{model_argu} #{bootstrap_argu} #{te_argu} #{add_argu}`
  end

  Thread.kill(thr) and puts if $? == 0

  `#{NW_TOPOLOGY} -Ib #{boottree_file} | ruby #{REORDER_NODE} -i - --ref #{ref_tree_file} > #{iqtree_outdir}/boot.bls`
  #`ruby #{REORDER_NODE} -i #{boottree_file} --ref #{ref_tree_file} > #{iqtree_outdir}/boot.bls`

  # for the ml tree (iqtree.treefile)
  `ruby #{REORDER_NODE} -i #{mltree_file} --ref #{ref_tree_file} > #{iqtree_outdir}/ml.bls`

  `#{RSCRIPT} #{FROM_BS_TO_HESSIAN} #{iqtree_outdir}/boot.bls #{iqtree_outdir}/hessian`

  mltree_file = File.join(iqtree_outdir, 'iqtree.treefile')
  create_inBV(mltree_file:mltree_file, mcmctree_outdir:mcmctree_outdir, inBV_file:inBV_file, iqtree_outdir:iqtree_outdir)
end


############################################################
inBV_files = get_files_under_folder(outdir, 'in.BV')

mcmctree_outdir = File.join(outdir, 'mcmctree')
mkdir_with_force(mcmctree_outdir)
`cat #{inBV_files} > #{mcmctree_outdir}/in.BV`

FileUtils.cp(ali_file, mcmctree_outdir)
`cp #{calib_tree_file} #{mcmctree_outdir}/species.trees`

prepare_paml_ctl(mcmctree_ctl_file, mcmctree_outdir, {'seqfile'=>'combined.phy', 'treefile'=>'species.trees'})

if is_run_mcmctree
  STDOUT.puts "Running MCMCTree ......"
  Dir.chdir(mcmctree_outdir)

  thr = processbar_for_bootstrapping(file:'mcmc.txt', b:get_no_iter(File.basename(mcmctree_ctl_file)))

  `#{MCMCTREE} #{File.basename(mcmctree_ctl_file)} > mcmctree.final`

  if $? == 0
    Thread.kill(thr) and puts
    `bash #{FIGTREE2NWK} -i FigTree.tre > figtree.nwk`
    puts "Done!" if $? == 0
  end
end