#! /bin/env ruby


###############################################################
require 'getoptlong'
require 'bio-nwk'


###############################################################
def get_ordered_tips(infile)
  ordered_tips = Array.new
  tree = nil

  in_fh = File.open(infile, 'r')
  in_fh.each_line do |line|
    line.chomp!
    line.scan(/([^(),: ]+):/) do |i|
      ordered_tips << i[0]
    end
    tree = getTreeObjFromNwkString(line)
    break
  end
  in_fh.close

  return([ordered_tips, tree])
end


def get_child(tree, node, bls, name2bl, name2order, order2names, is_name2order=false)
  if node != tree.root
    bls << tree.distance(node, tree.parent(node))
    tips = tree.tips(node)

    tip_names = tips.map{|i|i.name.gsub(' ', '_')}.sort
    complement_tip_names = (ORDERED_TIP_NAMES - tree.tips(node).map{|t|t.name.gsub(' ', '_')}).sort
    #p tree.tips(node).map{|t|t.name}

    name2bl[tip_names] = bls[-1]
    name2bl[complement_tip_names]= bls[-1]

    if is_name2order
      order = order2names.size # the 1st is zero
      name2order[tip_names] = order
      name2order[complement_tip_names] = order
      order2names[order] = [tip_names, complement_tip_names]
    end
  end

  if node.isTip?(tree)
    ;
  else
    children = tree.children(node)
    ordered_children = children.sort_by{|child| tree.tips(child).map{|i|ORDERED_TIP_NAMES.index(i.name.gsub(' ', '_'))}.min}
    #p [ordered_children, node]
    ordered_children.each do |i|
      get_child(tree, i, bls, name2bl, name2order, order2names, is_name2order)
    end
  end
end


###############################################################
infile = nil
ref_tree_file = nil

name2bl = Hash.new
name2order, order2names = [Hash.new, Hash.new]


###############################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
  ['--ref', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '-i'
      infile = value
    when '--ref'
      ref_tree_file = value
  end
end


###############################################################
ORDERED_TIP_NAMES, ref_tree = get_ordered_tips(ref_tree_file)

bls = Array.new
ref_tree = get_child(ref_tree, ref_tree.root, bls, name2bl, name2order, order2names, true)
NAME2ORDER = name2order
ORDER2NAMES = order2names


name2bl = Hash.new


###############################################################
trees = getTreeObjs(infile)

trees.each do |tree|
  bls = Array.new
  root = tree.root
  get_child(tree, root, bls, name2bl, name2order, order2names, false)

  bls2 = Array.new
  ORDER2NAMES.each_pair do |order, names|
    names = order2names[order]
    bls2 << names.map{|name| name2bl[name] }.compact[0]
    if names.map{|name| name2bl[name] }.any?{|i|i.nil?}
      p names
      p names.map{|name| name2bl[name] }
    end
    #p [bls[-1], names]
  end
  puts bls2.join(' ')
end


