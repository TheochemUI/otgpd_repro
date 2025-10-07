#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# --- Configuration ---
TAGS_DIR="apidocs/tags"
THEME_DIR="apidocs"

# External Tag File Configuration
TAG_URL="https://upload.cppreference.com/mwiki/images/f/f8/cppreference-doxygen-web.tag.xml"
TAG_FILE_NAME="cppreference-doxygen-web.tag.xml"
TAG_FILE_PATH="$TAGS_DIR/$TAG_FILE_NAME"

# Doxygen Theme Configuration
THEME_URL="https://github.com/HaoZeke/doxyYoda/releases/download/0.0.2/doxyYoda_0.0.2.tar.gz"
THEME_ARCHIVE_NAME="doxyYoda_0.0.2.tar.gz"
THEME_ARCHIVE_PATH="$THEME_DIR/$THEME_ARCHIVE_NAME"
THEME_EXTRACTED_DIR="doxyYoda_0.0.2"
THEME_EXTRACTED_PATH="$THEME_DIR/$THEME_EXTRACTED_DIR"

# --- Logic ---

# 1. Get External Doxygen Tags
echo "--- Checking for Doxygen tag file ---"
mkdir -p "$TAGS_DIR"

if [ -f "$TAG_FILE_PATH" ]; then
    echo "Tag file '$TAG_FILE_NAME' already exists. Skipping download."
else
    echo "Tag file not found. Downloading from $TAG_URL..."
    # Use curl with -L to follow redirects and -o to specify the output file
    curl -L "$TAG_URL" -o "$TAG_FILE_PATH"
    echo "Tag file downloaded successfully."
fi

echo "" # Add a newline for cleaner output

# 2. Get Doxygen Theme
echo "--- Checking for Doxygen theme ---"

if [ -d "$THEME_EXTRACTED_PATH" ]; then
    echo "Theme directory '$THEME_EXTRACTED_DIR' already exists. Skipping download and extraction."
else
    echo "Theme not found. Downloading and extracting..."
    # Use wget with -P to download the archive to the correct directory
    wget -P "$THEME_DIR" "$THEME_URL"

    # Extract the archive into the target directory
    tar -xf "$THEME_ARCHIVE_PATH" -C "$THEME_DIR"

    # Clean up the downloaded archive
    rm "$THEME_ARCHIVE_PATH"

    echo "Theme setup complete."
fi

echo ""
echo "Asset setup finished."
