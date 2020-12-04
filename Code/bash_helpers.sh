#!/bin/bash

alias build='docker build . -t bacpac_smartsims'

function rand-str {
    # Return random alpha-numeric string of given LENGTH
    #
    # Usage: VALUE=$(rand-str $LENGTH)
    #    or: VALUE=$(rand-str)

    local DEFAULT_LENGTH=64
    local LENGTH=${1:-$DEFAULT_LENGTH}

    tr -dc A-Za-z0-9 </dev/urandom | head -c $LENGTH
    # -dc: delete complementary set == delete all except given set
}

function box_out()
{
  local s=("$@") b w
  for l in "${s[@]}"; do
    ((w<${#l})) && { b="$l"; w="${#l}"; }
  done
  tput setaf 3
  echo " -${b//?/-}-
| ${b//?/ } |"
  for l in "${s[@]}"; do
    printf '| %s%*s%s |\n' "$(tput setaf 4)" "-$w" "$l" "$(tput setaf 3)"
  done
  echo "| ${b//?/ } |
 -${b//?/-}-"
  tput sgr 0
}

function run_rss(){
    export PW=`rand-str 16`
    
    box_out 'RStudio Running on Port 8787' 'User: rstudio' "Password: $PW" 'C-c C-c to kill container'

    docker run -v `pwd`/..:/home/rstudio/bacpac_smartsims \
       -p 8787:8787 \
       -e PASSWORD=$PW\
       -it bacpac_smartsims

}
