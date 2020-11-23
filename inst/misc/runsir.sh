#!/bin/sh
SIRBINPATH="$1"
JAVAPATH="$SIRBINPATH/../lib/runtime/bin/java"
LIBPATH="$SIRBINPATH/../lib"

shift

"$JAVAPATH" -Xms1G -XX:MaxRAMPercentage=85 --illegal-access=permit --add-opens java.base/java.lang=ALL-UNNAMED \
    --add-opens java.base/java.net=ALL-UNNAMED -Djava.library.path="$LIBPATH/native" \
    "-javaagent:$LIBPATH/app/agents-4.4.8.jar" -cp "$LIBPATH/app/*" \
    de.unijena.bioinf.ms.frontend.SiriusCLIApplication $@
