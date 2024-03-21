#! /bin/env ruby


def processbar(count, total)
  maxlen = 80
  prop = (count.to_f/total)*100
  jing = '#' * (prop*maxlen/100)
  printf "\r%-80s%s", jing, prop.round(2).to_s+'%'
end


