import os
import sys
import requests
import tempfile

def setup_ngl():
    try:
        # Get the project root directory (same as script location since script is in root)
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = script_dir
        ngl_assets_dir = os.path.join(project_root, "assets", "ngl_assets")
        
        print(f"Current working directory: {os.getcwd()}")
        print(f"Script directory: {script_dir}")
        print(f"Project root: {project_root}")
        print(f"Target directory: {ngl_assets_dir}")
        
        # Create the target directory with more verbose output
        try:
            os.makedirs(ngl_assets_dir, exist_ok=True)
            print(f"Directory created or verified: {ngl_assets_dir}")
            print(f"Directory exists: {os.path.exists(ngl_assets_dir)}")
            print(f"Directory is writable: {os.access(ngl_assets_dir, os.W_OK)}")
        except Exception as e:
            print(f"Error creating directory: {e}")
            return False

        ngl_url = "https://cdn.jsdelivr.net/npm/ngl@2.0.0-dev.37/dist/ngl.min.js"
        ngl_file_path = os.path.join(ngl_assets_dir, "ngl.min.js")
        
        # Try downloading to a temporary file first
        temp_file = None
        try:
            # Create a temporary file in the same directory to avoid cross-device issues
            temp_file = os.path.join(ngl_assets_dir, ".ngl.min.js.download")
            print(f"Downloading {ngl_url} to temporary file: {temp_file}")
            
            # Download the file
            response = requests.get(ngl_url, stream=True)
            response.raise_for_status()
            
            # Write to temporary file
            with open(temp_file, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:  # filter out keep-alive chunks
                        f.write(chunk)
            
            # Verify the temporary file
            if not os.path.exists(temp_file):
                print("Error: Temporary file was not created")
                return False
                
            temp_size = os.path.getsize(temp_file)
            print(f"Temporary file created: {temp_file} ({temp_size} bytes)")
            
            # Move the temporary file to the final location
            if os.path.exists(ngl_file_path):
                os.remove(ngl_file_path)
            os.rename(temp_file, ngl_file_path)
            
            # Final verification
            if os.path.exists(ngl_file_path):
                file_size = os.path.getsize(ngl_file_path)
                print(f"Successfully downloaded ngl.min.js to {ngl_file_path} ({file_size} bytes)")
                return True
            else:
                print("Error: File was not moved to final location")
                return False
                
        except Exception as e:
            print(f"Error during download: {e}")
            # Clean up temporary file if it exists
            if temp_file and os.path.exists(temp_file):
                try:
                    os.remove(temp_file)
                except:
                    pass
            return False
            
    except requests.exceptions.RequestException as e:
        print(f"Error downloading ngl.min.js: {e}")
        print("Please ensure you have an internet connection and try again.")
        return False
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = setup_ngl()
    sys.exit(0 if success else 1)
