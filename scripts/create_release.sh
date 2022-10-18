#!/bin/bash


release_name="0.3"
release_description="Dockerization"


release_notes="CALDER2 release: ${release_description} (${release_name})"


latest_tag=$(git describe --tags --abbrev=0)

if [[ "${release_name}" == "${latest_tag}" ]]
then
	echo "Overwriting exiting release (${release_name})"
	echo "---------------------------------------------"

	# Detelting Github release (should prompt a confirmation)
	gh release delete ${release_name}
	# Deleting local tag
	git tag -d ${release_name}
	# Pushing the removal of the tag to Github
	git push origin :${release_name}
fi

echo "Uploading new release (${release_name})"
echo "---------------------------------------------"
# Creating new tag
git tag -a ${release_name} -m ${release_description}
# Pushing new tag
git push origin ${release_name}
# Creating release from tag
gh release create ${release_name} --title "${release_description}" --notes "${release_notes}"