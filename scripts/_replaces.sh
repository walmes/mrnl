# Lista de pacotes que aparecem no conjunto de arquivos.
# ATTENTION: `(` é o literal e `\(` é o meta.
grep -h -E '(library|require)\(' *.R | sed 's/^.*(\(.*\))$/"\1",/g' | sort | uniq

# Show only the modified lines.
sed -n -e '/^\(require\|library\)/p' *.R
sed -n s/require/library/p script04-infograd.R

# Does the replacement in place.
sed --in-place s/require/library/g *.R
