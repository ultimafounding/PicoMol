import os
import sys
import requests

def setup_ngl():
    # Get the project root directory (one level up from the script location)
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ngl_assets_dir = os.path.join(project_root, "assets", "ngl_assets")
    
    # Create the ngl_assets directory if it doesn't exist
    os.makedirs(ngl_assets_dir, exist_ok=True)

    ngl_url = "https://cdn.jsdelivr.net/npm/ngl@2.0.0-dev.37/dist/ngl.min.js"
    ngl_file_path = os.path.join(ngl_assets_dir, "ngl.min.js")

    print(f"Downloading ngl.min.js from {ngl_url}...")
    try:
        response = requests.get(ngl_url, stream=True)
        response.raise_for_status()  # Raise an HTTPError for bad responses (4xx or 5xx)

        with open(ngl_file_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print(f"Successfully downloaded ngl.min.js to {ngl_assets_dir}")
        return True
    except requests.exceptions.RequestException as e:
        print(f"Error downloading ngl.min.js: {e}")
        print("Please ensure you have an internet connection and try again.")
        return False
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return False

if __name__ == "__main__":
    success = setup_ngl()
    sys.exit(0 if success else 1)
