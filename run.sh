#!/usr/bin/bash

echo "⚙️ Compiling..."
./src/compile.sh
echo "✅ Done!"

echo "🏃 Running..."
./data/main.out
echo "✅ Done!"