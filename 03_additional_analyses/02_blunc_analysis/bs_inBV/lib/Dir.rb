#! /usr/bin/env ruby

require 'find'


################################################################################
class Dir
    def self.mkdirs(path)
        if(!File.directory?(path))
            if(!mkdirs(File.dirname(path)))
                return false;
            end
            mkdir(path)
        end
        return true
    end
end


################################################################################
def mkdir_with_force(outdir, is_force=false, is_tolerate=false)
  if outdir.class != String
    raise "outdir wrong? Exiting ......"
  end

  if ! Dir.exists?(outdir)
    `mkdir -p #{outdir}`
  else
    if is_tolerate
      ;
    elsif is_force
      `rm -rf #{outdir}`
      `mkdir -p #{outdir}`
    else
      raise "The outdir #{outdir} has already existed!"
    end
  end
end


def read_infiles(indir, suffix='', is_all_subfolder=false)
  infiles = Array.new
  if ! is_all_subfolder
    Dir.foreach(indir) do |b|
      next if b =~ /^\./
      if suffix.is_a?(String)
        if suffix != ''
          next if b !~ /#{suffix}$/
        end
      elsif suffix.is_a?(Array)
        next unless suffix.any?{|i| b =~ /#{i}$/ }
      end
      infiles << File.join(indir, b)
    end
  else
    Find.find(indir) do |path|
      next if File.directory?(path)
      next if File.basename(path) =~ /^\./
      infiles << path if suffix.is_a?(String) ? path =~ /\.#{suffix}$/ : suffix.any?{|i| path =~ /#{i}$/ }
    end
  end
  return(infiles)
end


def getFilesBySuffices(indir, suffices)
  files = Array.new
  infiles = read_infiles(indir)
  infiles.each do |infile|
    if suffices.include?(File.extname(infile))
      files << infile
    end
  end
  return(files)
end


def get_file_path(file)
  path = File.symlink?(file) ? File.readlink(file) : file
  return(path)
end