FROM gcr.io/google-appengine/python

RUN apt-get update && apt-get install -y --no-install-recommends \
build-essential \
ca-certificates \
curl \
git \
libncursesw5-dev \
libncurses5-dev \
make \
zlib1g-dev \
libbz2-dev \
liblzma-dev \
fuse \
tabix \
graphviz \
libgraphviz-dev \
pkg-config \
&& rm -rf /var/lib/apt/lists/*

ENV BCFTOOLS_BIN="bcftools-1.4.tar.bz2" \
BCFTOOLS_PLUGINS="/usr/local/libexec/bcftools" \
BCFTOOLS_VERSION="1.4"

# Install BCFTools
RUN curl -fsSL https://github.com/samtools/bcftools/releases/download/$BCFTOOLS_VERSION/$BCFTOOLS_BIN -o /opt/$BCFTOOLS_BIN \
&& tar xvjf /opt/$BCFTOOLS_BIN -C /opt/ \
&& cd /opt/bcftools-$BCFTOOLS_VERSION \
&& make \
&& make install 2>&1

# Create a virtualenv for dependencies. This isolates these packages from
# system-level packages.
RUN virtualenv /env

# Setting these environment variables are the same as running
# source /env/bin/activate.
ENV VIRTUAL_ENV /env
ENV PATH /env/bin:$PATH

# Use python3.6!
RUN update-alternatives --install /usr/bin/python python /opt/python3.6/bin/python3.6 2 \
    && ln -f /opt/python3.6/bin/python3.6 /env/bin/python \
    && ln -f /opt/python3.6/bin/pip3.6 /env/bin/pip

# Copy the application's requirements.txt and run pip to install all
# dependencies into the virtualenv.
ADD requirements.txt /app/requirements.txt
RUN pip install -r /app/requirements.txt

# Add the application source code.
ADD . /app


CMD gunicorn -b :$PORT main:app
