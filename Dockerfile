# Start from a base image
FROM nextstrain/base:build-20240710T214955Z

# Install Python packages.
# Allow Snakemake to create subdirs in the user cache dir
# <https://github.com/nextstrain/ncov-ingest/pull/401>
RUN pip3 install \
    "altair[all]" \
    jupyterlab \
    notebook \
    seaborn \
    statsmodels \
 && rm -rf ~/.cache

# Install Java and C compiler.
RUN apt-get update && apt-get install -y --no-install-recommends \
    default-jre \
    clang && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install Julia.
# Recipe copied from official Julia Dockerfile:
# https://github.com/docker-library/julia/blob/master/1.10/bullseye/Dockerfile
ENV JULIA_PATH /usr/local/julia
ENV PATH $JULIA_PATH/bin:$PATH

# https://julialang.org/juliareleases.asc
# Julia (Binary signing key) <buildbot@julialang.org>
ENV JULIA_GPG 3673DF529D9049477F76B37566E3C7DC03D6E495

# https://julialang.org/downloads/
ENV JULIA_VERSION 1.10.4

RUN set -eux; \
	\
	savedAptMark="$(apt-mark showmanual)"; \
	apt-get update; \
	apt-get install -y --no-install-recommends \
		gnupg \
	; \
	rm -rf /var/lib/apt/lists/*; \
	\
# https://julialang.org/downloads/#julia-command-line-version
# https://julialang-s3.julialang.org/bin/checksums/julia-1.10.4.sha256
	arch="$(dpkg --print-architecture)"; \
	case "$arch" in \
		'amd64') \
			url='https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-1.10.4-linux-x86_64.tar.gz'; \
			sha256='079f61757c3b5b40d2ade052b3cc4816f50f7ef6df668825772562b3746adff1'; \
			;; \
		'arm64') \
			url='https://julialang-s3.julialang.org/bin/linux/aarch64/1.10/julia-1.10.4-linux-aarch64.tar.gz'; \
			sha256='ae4ae6ade84a103cdf30ce91c8d4035a0ef51c3e2e66f90a0c13abeb4e100fc4'; \
			;; \
		'i386') \
			url='https://julialang-s3.julialang.org/bin/linux/x86/1.10/julia-1.10.4-linux-i686.tar.gz'; \
			sha256='5771c6032b2be4442f4d4f4fc610eac09a48888ce1f4a41d10208e9d413f5054'; \
			;; \
		'ppc64el') \
			url='https://julialang-s3.julialang.org/bin/linux/ppc64le/1.10/julia-1.10.4-linux-ppc64le.tar.gz'; \
			sha256='0703f983894974491715e816a006e0f063966544023f470c94c71ef99dff9dba'; \
			;; \
		*) \
			echo >&2 "error: current architecture ($arch) does not have a corresponding Julia binary release"; \
			exit 1; \
			;; \
	esac; \
	\
	curl -fL -o julia.tar.gz.asc "$url.asc"; \
	curl -fL -o julia.tar.gz "$url"; \
	\
	echo "$sha256 *julia.tar.gz" | sha256sum --strict --check -; \
	\
	export GNUPGHOME="$(mktemp -d)"; \
	gpg --batch --keyserver keyserver.ubuntu.com --recv-keys "$JULIA_GPG"; \
	gpg --batch --verify julia.tar.gz.asc julia.tar.gz; \
	gpgconf --kill all; \
	rm -rf "$GNUPGHOME" julia.tar.gz.asc; \
	\
	mkdir "$JULIA_PATH"; \
	tar -xzf julia.tar.gz -C "$JULIA_PATH" --strip-components 1; \
	rm julia.tar.gz; \
	\
	apt-mark auto '.*' > /dev/null; \
	[ -z "$savedAptMark" ] || apt-mark manual $savedAptMark; \
	apt-get purge -y --auto-remove -o APT::AutoRemove::RecommendsImportant=false; \
	\
# smoke test
	julia --version

# Install TreeKnit.
RUN julia -e 'using Pkg; Pkg.add(["TreeKnit"]); Pkg.build("TreeKnit");'
