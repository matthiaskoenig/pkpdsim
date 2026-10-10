"""The creators and the terms of use of a model.

A model sets them in its `sbmlutils.factory.Model`: the creators as
`creators`, the terms of use appended to its `notes`, behind the description
of the model. A model with further creators extends the list, e.g.
`creators=[*creators, Creator(...)]`.
"""

from datetime import datetime

from sbmlutils.factory import Creator

creators = [
    Creator(
        familyName="König",
        givenName="Matthias",
        email="koenigmx@hu-berlin.de",
        organization="Humboldt-University Berlin, Institute for Theoretical Biology",
        site="https://livermetabolism.com",
    ),
]
"""The creators of the models, the first one is named in the terms of use."""

_TERMS_OF_USE = """
    The content of this model has been carefully created in a manual research effort.

    ## Terms of use
    Copyright © {year} <a href="{site}" title="{given_name} {family_name}" target="_blank">{given_name} {family_name}</a>.

    <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">
    <img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

    Redistribution and use of any part of this model, with or without modification, are permitted provided
    that the following conditions are met:

    - Redistributions of this SBML file must retain the above copyright notice, this list of conditions and the following disclaimer.
    - Redistributions in a different form must reproduce the above copyright notice, this list of conditions
      and the following disclaimer in the documentation and/or other materials provided
      with the distribution.

    This model is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
    implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
"""
"""The template of the terms of use, `terms_of_use` fills in the year and the creator."""

terms_of_use = _TERMS_OF_USE.format(
    year=datetime.now().year,
    given_name=creators[0].givenName,
    family_name=creators[0].familyName,
    site=creators[0].site,
)
"""The terms of use of a model in markdown, with the copyright of the current year."""
